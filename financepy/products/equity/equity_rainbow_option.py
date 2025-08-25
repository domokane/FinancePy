##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from typing import List
from enum import Enum

from math import exp, log, sqrt
import numpy as np

from ...utils.math import normcdf
from ...utils.math import M
from ...utils.global_vars import G_DAYS_IN_YEARS
from ...utils.error import FinError
from ...models.gbm_process_simulator import get_assets_paths
from ...products.equity.equity_option import EquityOption
from ...market.curves.discount_curve import DiscountCurve
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.date import Date


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
    stock_prices,
    discount_curve,
    dividend_curves,
    volatilities,
    betas,
    num_assets,
    payoff_type,
    payoff_params,
    num_paths=10000,
    seed=4242,
):

    np.random.seed(seed)

    df = discount_curve.df_t(t)
    r = -log(df) / t

    qs = []
    for curve in dividend_curves:
        dq = curve.df_t(t)
        q = -np.log(dq) / t
        qs.append(q)

    qs = np.array(qs)

    mus = r - qs

    _, s_all = get_assets_paths(
        num_assets,
        num_paths,
        t,
        mus,
        stock_prices,
        volatilities,
        betas,
        seed,
    )

    payoff = payoff_value(s_all, payoff_type.value, payoff_params)
    payoff = np.mean(payoff)
    v = payoff * exp(-r * t)
    return v


########################################################################################


class EquityRainbowOption(EquityOption):

    def __init__(
        self,
        expiry_dt: Date,
        payoff_type: EquityRainbowOptionTypes,
        payoff_params: List[float],
        num_assets: int,
    ):

        check_argument_types(self.__init__, locals())

        self._validate_payoff(payoff_type, payoff_params, num_assets)

        self.expiry_dt = expiry_dt
        self.payoff_type = payoff_type
        self.payoff_params = payoff_params
        self.num_assets = num_assets

    ###########################################################################

    def _validate(self, stock_prices, dividend_curves, volatilities, betas):

        if len(stock_prices) != self.num_assets:
            raise FinError(
                "Stock prices must be a vector of length " + str(self.num_assets)
            )

        if len(dividend_curves) != self.num_assets:
            raise FinError(
                "Dividend discount must be a vector of length " + str(self.num_assets)
            )

        if len(volatilities) != self.num_assets:
            raise FinError(
                "Volatilities must be a vector of length " + str(self.num_assets)
            )

        if len(betas) != self.num_assets:
            raise FinError("Betas must be a vector of length " + str(self.num_assets))

    ###########################################################################

    def _validate_payoff(self, payoff_type, payoff_params, num_assets):

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
            raise FinError(
                "Number of parameters required for "
                + str(payoff_type)
                + " must be "
                + str(num_params)
            )

        if (
            payoff_type == EquityRainbowOptionTypes.CALL_ON_NTH
            or payoff_type == EquityRainbowOptionTypes.PUT_ON_NTH
        ):
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

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self.expiry_dt:
            raise FinError("Valuation date after expiry date.")

        if discount_curve.value_dt != value_dt:
            raise FinError(
                "Discount Curve valuation date not same as option value date"
            )

        if self.num_assets != 2:
            raise FinError("Analytical results for two assets only.")

        if corr_matrix.ndim != 2:
            raise FinError("Corr matrix must be of size 2x2")

        if corr_matrix.shape[0] != 2:
            raise FinError("Corr matrix must be of size 2x2")

        if corr_matrix.shape[1] != 2:
            raise FinError("Corr matrix must be of size 2x2")

        if value_dt > self.expiry_dt:
            raise FinError("Value date after expiry date.")

        # Use result by Stulz (1982) given by Haug Page 211
        t = (self.expiry_dt - value_dt) / G_DAYS_IN_YEARS
        r = discount_curve.zero_rate(self.expiry_dt)

        q1 = dividend_curves[0].zero_rate(self.expiry_dt)
        q2 = dividend_curves[1].zero_rate(self.expiry_dt)

        dividend_yields = [q1, q2]

        self._validate(stock_prices, dividend_yields, volatilities, corr_matrix)

        #        q1 = dividend_yields[0]
        #        q2 = dividend_yields[1]

        rho = corr_matrix[0][1]
        s1 = stock_prices[0]
        s2 = stock_prices[1]
        b1 = r - q1
        b2 = r - q2
        v1 = volatilities[0]
        v2 = volatilities[1]
        k = self.payoff_params[0]

        v = sqrt(v1 * v1 + v2 * v2 - 2 * rho * v1 * v2)
        d = (log(s1 / s2) + (b1 - b2 + v * v / 2) * t) / v / sqrt(t)
        y1 = (log(s1 / k) + (b1 + v1 * v1 / 2) * t) / v1 / sqrt(t)
        y2 = (log(s2 / k) + (b2 + v2 * v2 / 2) * t) / v2 / sqrt(t)
        rho1 = (v1 - rho * v2) / v
        rho2 = (v2 - rho * v1) / v
        dq1 = exp(-q1 * t)
        dq2 = exp(-q2 * t)
        df = exp(-r * t)

        if self.payoff_type == EquityRainbowOptionTypes.CALL_ON_MAXIMUM:
            v = (
                s1 * dq1 * M(y1, d, rho1)
                + s2 * dq2 * M(y2, -d + v * sqrt(t), rho2)
                - k * df * (1.0 - M(-y1 + v1 * sqrt(t), -y2 + v2 * sqrt(t), rho))
            )
        elif self.payoff_type == EquityRainbowOptionTypes.CALL_ON_MINIMUM:
            v = (
                s1 * dq1 * M(y1, -d, -rho1)
                + s2 * dq2 * M(y2, d - v * sqrt(t), -rho2)
                - k * df * M(y1 - v1 * sqrt(t), y2 - v2 * sqrt(t), rho)
            )
        elif self.payoff_type == EquityRainbowOptionTypes.PUT_ON_MAXIMUM:
            cmax1 = (
                s2 * dq2 + s1 * dq1 * normcdf(d) - s2 * dq2 * normcdf(d - v * sqrt(t))
            )
            cmax2 = (
                s1 * dq1 * M(y1, d, rho1)
                + s2 * dq2 * M(y2, -d + v * sqrt(t), rho2)
                - k * df * (1.0 - M(-y1 + v1 * sqrt(t), -y2 + v2 * sqrt(t), rho))
            )
            v = k * df - cmax1 + cmax2
        elif self.payoff_type == EquityRainbowOptionTypes.PUT_ON_MINIMUM:
            cmin1 = (
                s1 * dq1 - s1 * dq1 * normcdf(d) + s2 * dq2 * normcdf(d - v * sqrt(t))
            )
            cmin2 = (
                s1 * dq1 * M(y1, -d, -rho1)
                + s2 * dq2 * M(y2, d - v * sqrt(t), -rho2)
                - k * df * M(y1 - v1 * sqrt(t), y2 - v2 * sqrt(t), rho)
            )
            v = k * df - cmin1 + cmin2
        else:
            raise FinError("Unsupported Rainbow option type")

        return v

    ###########################################################################

    def value_mc(
        self,
        value_dt,
        stock_prices,
        discount_curve,
        dividend_curves,
        volatilities,
        corr_matrix,
        num_paths=10000,
        seed=4242,
    ):

        self._validate(stock_prices, dividend_curves, volatilities, corr_matrix)

        if value_dt > self.expiry_dt:
            raise FinError("Value date after expiry date.")

        t = (self.expiry_dt - value_dt) / G_DAYS_IN_YEARS

        v = value_mc_fast(
            t,
            stock_prices,
            discount_curve,
            dividend_curves,
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

        s = label_to_string("OBJECT TYPE", type(self).__name__)
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
