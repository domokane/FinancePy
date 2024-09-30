##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO: Consider risk management
# TODO: Consider using Sobol
# TODO: Consider allowing weights on the individual basket assets
# TODO: Extend monte carlo to handle American options

import numpy as np

from ...utils.global_vars import g_days_in_year

from ...models.gbm_process_simulator import get_assets_paths

from ...utils.error import FinError
from ...utils.global_types import OptionTypes
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.helpers import _func_name
from ...utils.date import Date
from ...market.curves.discount_curve import DiscountCurve

from ...utils.math import N


###############################################################################


class EquityBasketOption:
    """A EquityBasketOption is a contract to buy a put or a call option on
    an equally weighted portfolio of different stocks, each with its own price,
    volatility and dividend yield. An analytical and monte-carlo pricing model
    have been implemented for a European style option."""

    def __init__(
        self,
        expiry_dt: Date,
        strike_price: float,
        option_type: OptionTypes,
        num_assets: int,
    ):
        """Define the EquityBasket option by specifying its expiry date,
        its strike price, whether it is a put or call, and the number of
        underlying stocks in the basket."""

        check_argument_types(self.__init__, locals())

        self.expiry_dt = expiry_dt
        self.strike_price = float(strike_price)
        self.option_type = option_type
        self.num_assets = num_assets

    ###########################################################################

    def _validate(
        self, stock_prices, dividend_yields, volatilities, correlations
    ):

        if len(stock_prices) != self.num_assets:
            raise FinError(
                "Stock prices must have a length " + str(self.num_assets)
            )

        if len(dividend_yields) != self.num_assets:
            raise FinError(
                "Dividend yields must have a length " + str(self.num_assets)
            )

        if len(volatilities) != self.num_assets:
            raise FinError(
                "Volatilities must have a length " + str(self.num_assets)
            )

        if correlations.ndim != 2:
            raise FinError("Correlation must be a 2D matrix ")

        if correlations.shape[0] != self.num_assets:
            raise FinError(
                "Correlation cols must have a length " + str(self.num_assets)
            )

        if correlations.shape[1] != self.num_assets:
            raise FinError(
                "correlation rows must have a length " + str(self.num_assets)
            )

        for i in range(0, self.num_assets):
            if correlations[i, i] != 1.0:
                raise FinError("Corr matrix must have 1.0 on the diagonal")

            for j in range(0, i):
                if abs(correlations[i, j]) > 1.0:
                    raise FinError("Correlations must be [-1, +1]")

                if abs(correlations[j, i]) > 1.0:
                    raise FinError("Correlations must be [-1, +1]")

                if correlations[i, j] != correlations[j, i]:
                    raise FinError("Correlation matrix must be symmetric")

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        stock_prices: np.ndarray,
        discount_curve: DiscountCurve,
        dividend_curves: list,
        volatilities: np.ndarray,
        correlations: np.ndarray,
    ):
        """Basket valuation using a moment matching method to approximate the
        effective variance of the underlying basket value. This approach is
        able to handle a full rank correlation structure between the individual
        assets."""

        t_exp = (self.expiry_dt - value_dt) / g_days_in_year

        if value_dt > self.expiry_dt:
            raise FinError("Value date after expiry date.")

        qs = []
        for curve in dividend_curves:
            q = curve.cc_rate(self.expiry_dt)
            qs.append(q)

        v = volatilities
        s = stock_prices

        self._validate(stock_prices, qs, volatilities, correlations)

        a = np.ones(self.num_assets) * (1.0 / self.num_assets)

        r = discount_curve.cc_rate(self.expiry_dt)

        smean = 0.0
        for ia in range(0, self.num_assets):
            smean = smean + s[ia] * a[ia]

        ln_s0_k = np.log(smean / self.strike_price)
        sqrt_t = np.sqrt(t_exp)

        # Moment matching - starting with dividend
        qnum = 0.0
        qden = 0.0
        for ia in range(0, self.num_assets):
            qnum = qnum + a[ia] * s[ia] * np.exp(-qs[ia] * t_exp)
            qden = qden + a[ia] * s[ia]
        qhat = -np.log(qnum / qden) / t_exp

        # Moment matching - matching volatility
        vnum = 0.0
        for ia in range(0, self.num_assets):
            for ja in range(0, ia):
                rho_sigma_sigma = v[ia] * v[ja] * correlations[ia, ja]
                exp_term = (qs[ia] + qs[ja] - rho_sigma_sigma) * t_exp
                vnum = vnum + a[ia] * a[ja] * s[ia] * s[ja] * np.exp(-exp_term)

        vnum *= 2.0

        for ia in range(0, self.num_assets):
            rho_sigma_sigma = v[ia] ** 2
            exp_term = (2.0 * qs[ia] - rho_sigma_sigma) * t_exp
            vnum = vnum + ((a[ia] * s[ia]) ** 2) * np.exp(-exp_term)

        vhat2 = np.log(vnum / qnum / qnum) / t_exp

        den = np.sqrt(vhat2) * sqrt_t
        mu = r - qhat
        d1 = (ln_s0_k + (mu + vhat2 / 2.0) * t_exp) / den
        d2 = (ln_s0_k + (mu - vhat2 / 2.0) * t_exp) / den

        if self.option_type == OptionTypes.EUROPEAN_CALL:
            v = smean * np.exp(-qhat * t_exp) * N(d1)
            v = v - self.strike_price * np.exp(-r * t_exp) * N(d2)
        elif self.option_type == OptionTypes.EUROPEAN_PUT:
            v = self.strike_price * np.exp(-r * t_exp) * N(-d2)
            v = v - smean * np.exp(-qhat * t_exp) * N(-d1)
        else:
            raise FinError("Unknown option type")

        return v

    ###########################################################################

    def value_mc(
        self,
        value_dt: Date,
        stock_prices: np.ndarray,
        discount_curve: DiscountCurve,
        dividend_curves: list,
        volatilities: np.ndarray,
        corr_matrix: np.ndarray,
        num_paths: int = 10000,
        seed: int = 4242,
    ):
        """Valuation of the EquityBasketOption using a Monte-Carlo simulation
        of stock prices assuming a GBM distribution. Cholesky decomposition is
        used to handle a full rank correlation structure between the individual
        assets. The num_paths and seed are pre-set to default values but can be
        overwritten."""

        check_argument_types(getattr(self, _func_name(), None), locals())

        if value_dt > self.expiry_dt:
            raise FinError("Value date after expiry date.")

        t_exp = (self.expiry_dt - value_dt) / g_days_in_year

        dividend_yields = []
        for curve in dividend_curves:
            dq = curve.df(self.expiry_dt)
            q = -np.log(dq) / t_exp
            dividend_yields.append(q)

        self._validate(
            stock_prices, dividend_yields, volatilities, corr_matrix
        )

        num_assets = len(stock_prices)

        df = discount_curve.df(self.expiry_dt)
        r = -np.log(df) / t_exp

        mus = r - dividend_yields
        k = self.strike_price

        np.random.seed(seed)

        t_all, s_all = get_assets_paths(
            num_assets,
            num_paths,
            t_exp,
            mus,
            stock_prices,
            volatilities,
            corr_matrix,
            seed,
        )

        if self.option_type == OptionTypes.EUROPEAN_CALL:
            payoff = np.maximum(np.mean(s_all, axis=0) - k, 0.0)
        elif self.option_type == OptionTypes.EUROPEAN_PUT:
            payoff = np.maximum(k - np.mean(s_all, axis=0), 0.0)
        else:
            raise FinError("Unknown option type.")

        payoff = np.mean(payoff)
        v = payoff * np.exp(-r * t_exp)
        return v

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("STRIKE PRICE", self.strike_price)
        s += label_to_string("OPTION TYPE", self.option_type)
        s += label_to_string("NUM ASSETS", self.num_assets, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


###############################################################################
