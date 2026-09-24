##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from typing import List
from enum import Enum

from math import exp, log, sqrt
import numpy as np

from ...utils.math import normcdf
from ...utils.math import M
from ...utils.error import FinError
from ...models.gbm_process_simulator import get_assets_paths
from ...products.equity.equity_option import EquityOption
from ...market.curves.discount_curve import DiscountCurve
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.date import Date

from ...utils.check_values import check_curve_dt
from ...utils.check_values import check_corr_matrix
from ...utils.check_values import check_volatility
from ...utils.check_values import check_stock_price
from ...utils.helpers import option_years


class EquityRainbowOptionTypes(Enum):
    CALL_ON_MAXIMUM = 1
    PUT_ON_MAXIMUM = 2
    CALL_ON_MINIMUM = 3
    PUT_ON_MINIMUM = 4
    CALL_ON_NTH = 5  # MAX(NTH(S1,S2,...,SN)-K,0)
    PUT_ON_NTH = 6  # MAX(K-NTH(S1,S2,...,SN),0)


########################################################################################


def payoff_value(s, payoff_type_value, payoff_params):

    if payoff_type_value == EquityRainbowOptionTypes.CALL_ON_MINIMUM.value:
        k = payoff_params[0]
        # average on asset
        payoff = np.maximum(np.min(s, axis=0) - k, 0.0)
    elif payoff_type_value == EquityRainbowOptionTypes.CALL_ON_MAXIMUM.value:
        k = payoff_params[0]
        # average on asset
        payoff = np.maximum(np.max(s, axis=0) - k, 0.0)
    elif payoff_type_value == EquityRainbowOptionTypes.PUT_ON_MINIMUM.value:
        k = payoff_params[0]
        # average on asset
        payoff = np.maximum(k - np.min(s, axis=0), 0.0)
    elif payoff_type_value == EquityRainbowOptionTypes.PUT_ON_MAXIMUM.value:
        k = payoff_params[0]
        # average on asset
        payoff = np.maximum(k - np.max(s, axis=0), 0.0)
    elif payoff_type_value == EquityRainbowOptionTypes.CALL_ON_NTH.value:
        n = payoff_params[0]
        k = payoff_params[1]
        # sort on asset
        ssorted = np.sort(s, axis=0)
        assetn = ssorted[-n, :]
        payoff = np.maximum(assetn - k, 0.0)
    elif payoff_type_value == EquityRainbowOptionTypes.PUT_ON_NTH.value:
        n = payoff_params[0]
        k = payoff_params[1]
        # sort on asset
        ssorted = np.sort(s, axis=0)
        assetn = ssorted[-n, :]
        payoff = np.maximum(k - assetn, 0.0)
    else:
        raise FinError("Unknown payoff type")

    return payoff


########################################################################################


def value_mc_fast(
    t,
    stock_prices: float,
    r,
    qs,
    volatilities: np.ndarray,
    corr_matrix: np.ndarray,
    num_assets: int,
    payoff_type,
    payoff_params,
    num_paths: int,
    seed=4242,
):

    mus = r - qs

    _, s_all = get_assets_paths(
        num_assets,
        num_paths,
        t,
        mus,
        stock_prices,
        volatilities,
        corr_matrix,
        seed,
    )

    payoff = payoff_value(s_all, payoff_type.value, payoff_params)
    payoff = np.mean(payoff)
    v = payoff * exp(-r * t)
    return v

########################################################################################


def value_mc_fast_cv(
    t,
    stock_prices: float,
    r,
    qs,
    volatilities: np.ndarray,
    corr_matrix: np.ndarray,
    num_assets: int,
    payoff_type,
    payoff_params,
    num_paths: int,
    seed=4242,
):
    """Monte Carlo rainbow option valuation using antithetic paths and
    individual vanilla-option control variates."""

    if num_assets != 2:
        raise FinError(
            "Control variate currently implemented for two assets only."
        )

    mus = r - qs

    _, s_all = get_assets_paths(
        num_assets,
        num_paths,
        t,
        mus,
        stock_prices,
        volatilities,
        corr_matrix,
        seed,
    )

    # Rainbow payoff.
    x = payoff_value(
        s_all,
        payoff_type.value,
        payoff_params,
    )

    k = payoff_params[0]
    df = exp(-r * t)

    # --------------------------------------------------------------
    # Vanilla control payoffs
    # --------------------------------------------------------------

    if payoff_type in (
        EquityRainbowOptionTypes.CALL_ON_MAXIMUM,
        EquityRainbowOptionTypes.CALL_ON_MINIMUM,
    ):

        y1 = np.maximum(s_all[0, :] - k, 0.0)
        y2 = np.maximum(s_all[1, :] - k, 0.0)

    elif payoff_type in (
        EquityRainbowOptionTypes.PUT_ON_MAXIMUM,
        EquityRainbowOptionTypes.PUT_ON_MINIMUM,
    ):

        y1 = np.maximum(k - s_all[0, :], 0.0)
        y2 = np.maximum(k - s_all[1, :], 0.0)

    else:
        raise FinError(
            "Control variate not supported for this payoff type."
        )

    # Discount all simulated payoffs.
    x = df * x
    y1 = df * y1
    y2 = df * y2

    # --------------------------------------------------------------
    # Exact expectations of the vanilla controls
    #
    # These are ordinary Black-Scholes values.
    # --------------------------------------------------------------

    s1 = stock_prices[0]
    s2 = stock_prices[1]

    q1 = qs[0]
    q2 = qs[1]

    v1 = volatilities[0]
    v2 = volatilities[1]

    sqrt_t = sqrt(t)

    d11 = (
        log(s1 / k)
        + (r - q1 + 0.5 * v1 * v1) * t
    ) / (v1 * sqrt_t)

    d12 = d11 - v1 * sqrt_t

    d21 = (
        log(s2 / k)
        + (r - q2 + 0.5 * v2 * v2) * t
    ) / (v2 * sqrt_t)

    d22 = d21 - v2 * sqrt_t

    dq1 = exp(-q1 * t)
    dq2 = exp(-q2 * t)

    if payoff_type in (
        EquityRainbowOptionTypes.CALL_ON_MAXIMUM,
        EquityRainbowOptionTypes.CALL_ON_MINIMUM,
    ):

        exact_y1 = (
            s1 * dq1 * normcdf(d11)
            - k * df * normcdf(d12)
        )

        exact_y2 = (
            s2 * dq2 * normcdf(d21)
            - k * df * normcdf(d22)
        )

    else:

        exact_y1 = (
            k * df * normcdf(-d12)
            - s1 * dq1 * normcdf(-d11)
        )

        exact_y2 = (
            k * df * normcdf(-d22)
            - s2 * dq2 * normcdf(-d21)
        )

    # --------------------------------------------------------------
    # Antithetic pair averages
    #
    # get_assets_paths stores:
    #
    #     Z, -Z, Z, -Z, ...
    #
    # so each adjacent pair is one independent observation.
    # --------------------------------------------------------------

    x_pair = 0.5 * (
        x[0::2] + x[1::2]
    )

    y1_pair = 0.5 * (
        y1[0::2] + y1[1::2]
    )

    y2_pair = 0.5 * (
        y2[0::2] + y2[1::2]
    )

    # --------------------------------------------------------------
    # Estimate optimal control-variate coefficients
    #
    # beta = Cov(Y,Y)^(-1) Cov(Y,X)
    # --------------------------------------------------------------

    x_mean = np.mean(x_pair)

    y1_mean = np.mean(y1_pair)
    y2_mean = np.mean(y2_pair)

    dx = x_pair - x_mean
    dy1 = y1_pair - y1_mean
    dy2 = y2_pair - y2_mean

    var_y1 = np.mean(dy1 * dy1)
    var_y2 = np.mean(dy2 * dy2)
    cov_y1_y2 = np.mean(dy1 * dy2)

    cov_y1_x = np.mean(dy1 * dx)
    cov_y2_x = np.mean(dy2 * dx)

    cov_yy = np.array(
        [
            [var_y1, cov_y1_y2],
            [cov_y1_y2, var_y2],
        ]
    )

    cov_yx = np.array(
        [
            cov_y1_x,
            cov_y2_x,
        ]
    )

    beta = np.linalg.solve(
        cov_yy,
        cov_yx,
    )

    # --------------------------------------------------------------
    # Control-variate estimator
    # --------------------------------------------------------------

    x_cv = (
        x_pair
        - beta[0] * (y1_pair - exact_y1)
        - beta[1] * (y2_pair - exact_y2)
    )

    return np.mean(x_cv)

########################################################################################


class EquityRainbowOption(EquityOption):

    def __init__(
        self,
        expiry_dt: Date,
        payoff_type: EquityRainbowOptionTypes,
        payoff_params: List[float],
        num_assets: int,
    ) -> None:

        check_argument_types(self.__init__, locals())

        self._validate_payoff(payoff_type, payoff_params, num_assets)

        self.expiry_dt = expiry_dt
        self.payoff_type = payoff_type
        self.payoff_params = payoff_params
        self.num_assets = num_assets

    ###########################################################################

    def _validate_payoff(self, payoff_type, payoff_params, num_assets: int):

        num_params = 0

        if payoff_type == EquityRainbowOptionTypes.CALL_ON_MINIMUM:
            num_params = 1
        elif payoff_type == EquityRainbowOptionTypes.CALL_ON_MAXIMUM:
            num_params = 1
        elif payoff_type == EquityRainbowOptionTypes.PUT_ON_MINIMUM:
            num_params = 1
        elif payoff_type == EquityRainbowOptionTypes.PUT_ON_MAXIMUM:
            num_params = 1
        elif payoff_type == EquityRainbowOptionTypes.CALL_ON_NTH:
            num_params = 2
        elif payoff_type == EquityRainbowOptionTypes.PUT_ON_NTH:
            num_params = 2
        else:
            raise FinError("Unknown payoff type")

        if len(payoff_params) != num_params:
            raise FinError("Number of parameters required for " + str(payoff_type) + " must be " + str(num_params))

        if payoff_type == EquityRainbowOptionTypes.CALL_ON_NTH or payoff_type == EquityRainbowOptionTypes.PUT_ON_NTH:
            n = payoff_params[0]
            if n < 1 or n > num_assets:
                raise FinError("Nth parameter must be 1 to " + str(num_assets))

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        stock_prices: np.ndarray,
        discount_curve: DiscountCurve,
        dividend_curves: list,
        volatilities: np.ndarray,
        corr_matrix: np.ndarray,
    ):

        if self.num_assets != 2:
            raise FinError("Analytical results for two assets only.")

        t_exp = option_years(value_dt, self.expiry_dt)

        check_stock_price(stock_prices)
        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, *dividend_curves)
        check_corr_matrix(corr_matrix, self.num_assets)
        check_volatility(volatilities)

        # Use result by Stulz (1982) given by Haug Page 211

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q1 = dividend_curves[0].zero_rate_cc(self.expiry_dt)
        q2 = dividend_curves[1].zero_rate_cc(self.expiry_dt)

        rho = corr_matrix[0][1]
        s1 = stock_prices[0]
        s2 = stock_prices[1]
        b1 = r - q1
        b2 = r - q2
        v1 = volatilities[0]
        v2 = volatilities[1]
        k = self.payoff_params[0]

        v_sq = (v1 * v1 + v2 * v2 - 2.0 * rho * v1 * v2)

        if v_sq <= 1e-14:
            raise FinError(
                "Rainbow analytic formula is singular "
                "for this correlation/volatility combination."
            )

        v = sqrt(v_sq)

        d = (log(s1 / s2) + (b1 - b2 + v * v / 2) * t_exp) / v / sqrt(t_exp)
        y1 = (log(s1 / k) + (b1 + v1 * v1 / 2) * t_exp) / v1 / sqrt(t_exp)
        y2 = (log(s2 / k) + (b2 + v2 * v2 / 2) * t_exp) / v2 / sqrt(t_exp)
        rho1 = (v1 - rho * v2) / v
        rho2 = (v2 - rho * v1) / v
        dq1 = exp(-q1 * t_exp)
        dq2 = exp(-q2 * t_exp)
        df = exp(-r * t_exp)

        if self.payoff_type == EquityRainbowOptionTypes.CALL_ON_MAXIMUM:
            v = (
                s1 * dq1 * M(y1, d, rho1)
                + s2 * dq2 * M(y2, -d + v * sqrt(t_exp), rho2)
                - k * df * (1.0 - M(-y1 + v1 * sqrt(t_exp), -y2 + v2 * sqrt(t_exp), rho))
            )
        elif self.payoff_type == EquityRainbowOptionTypes.CALL_ON_MINIMUM:
            v = (
                s1 * dq1 * M(y1, -d, -rho1)
                + s2 * dq2 * M(y2, d - v * sqrt(t_exp), -rho2)
                - k * df * M(y1 - v1 * sqrt(t_exp), y2 - v2 * sqrt(t_exp), rho)
            )
        elif self.payoff_type == EquityRainbowOptionTypes.PUT_ON_MAXIMUM:
            cmax1 = s2 * dq2 + s1 * dq1 * normcdf(d) - s2 * dq2 * normcdf(d - v * sqrt(t_exp))
            cmax2 = (
                s1 * dq1 * M(y1, d, rho1)
                + s2 * dq2 * M(y2, -d + v * sqrt(t_exp), rho2)
                - k * df * (1.0 - M(-y1 + v1 * sqrt(t_exp), -y2 + v2 * sqrt(t_exp), rho))
            )
            v = k * df - cmax1 + cmax2
        elif self.payoff_type == EquityRainbowOptionTypes.PUT_ON_MINIMUM:
            cmin1 = s1 * dq1 - s1 * dq1 * normcdf(d) + s2 * dq2 * normcdf(d - v * sqrt(t_exp))
            cmin2 = (
                s1 * dq1 * M(y1, -d, -rho1)
                + s2 * dq2 * M(y2, d - v * sqrt(t_exp), -rho2)
                - k * df * M(y1 - v1 * sqrt(t_exp), y2 - v2 * sqrt(t_exp), rho)
            )
            v = k * df - cmin1 + cmin2
        else:
            raise FinError("Unsupported Rainbow option type")

        return v

    ###########################################################################

    def value_mc(
        self,
        value_dt: Date,
        stock_prices: float,
        discount_curve: DiscountCurve,
        dividend_curves: list[DiscountCurve],
        volatilities: np.ndarray,
        corr_matrix: np.ndarray,
        num_paths: int,
        seed=4242,
    ):

        t_exp = option_years(value_dt, self.expiry_dt)
        check_stock_price(stock_prices)
        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, *dividend_curves)
        check_corr_matrix(corr_matrix, self.num_assets)
        check_volatility(volatilities)

        r = discount_curve.zero_rate_cc(self.expiry_dt)

        qs = []
        for curve in dividend_curves:
            q = curve.zero_rate_cc(self.expiry_dt)
            qs.append(q)

        qs = np.array(qs)

        v = value_mc_fast(
            t_exp,
            stock_prices,
            r,
            qs,
            volatilities,
            corr_matrix,
            self.num_assets,
            self.payoff_type,
            self.payoff_params,
            num_paths,
            seed,
        )

        return v

    ###########################################################################

    def value_mc_cv(
        self,
        value_dt: Date,
        stock_prices: float,
        discount_curve: DiscountCurve,
        dividend_curves: list[DiscountCurve],
        volatilities: np.ndarray,
        corr_matrix: np.ndarray,
        num_paths: int,
        seed=4242,
    ):

        t_exp = option_years(
            value_dt,
            self.expiry_dt,
        )

        check_stock_price(stock_prices)
        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, *dividend_curves)
        check_corr_matrix(corr_matrix, self.num_assets)
        check_volatility(volatilities)

        r = discount_curve.zero_rate_cc(
            self.expiry_dt
        )

        qs = []

        for curve in dividend_curves:
            q = curve.zero_rate_cc(
                self.expiry_dt
            )
            qs.append(q)

        qs = np.array(qs)

        v = value_mc_fast_cv(
            t_exp,
            stock_prices,
            r,
            qs,
            volatilities,
            corr_matrix,
            self.num_assets,
            self.payoff_type,
            self.payoff_params,
            num_paths,
            seed,
        )

        return v

    ###########################################################################

    def __repr__(self):

        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("PAYOFF TYPE", self.payoff_type)
        s += label_to_string("PAYOFF PARAMS", self.payoff_params)
        s += label_to_string("NUM ASSETS TYPE", self.num_assets, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
