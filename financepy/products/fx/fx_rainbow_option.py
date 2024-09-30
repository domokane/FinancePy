##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from enum import Enum

import numpy as np
from typing import List

from ...utils.date import Date
from ...utils.math import N, M
from ...utils.global_vars import g_days_in_year
from ...utils.error import FinError
from ...models.gbm_process_simulator import get_assets_paths_times
from ...products.fx.fx_option import FXOption

from ...utils.helpers import check_argument_types


###############################################################################


class FXRainbowOptionTypes(Enum):
    CALL_ON_MAXIMUM = 1
    PUT_ON_MAXIMUM = 2
    CALL_ON_MINIMUM = 3
    PUT_ON_MINIMUM = 4
    CALL_ON_NTH = 5  # MAX(NTH(S1,S2,...,SN)-K,0)
    PUT_ON_NTH = 6  # MAX(K-NTH(S1,S2,...,SN),0)


###############################################################################


def payoff_value(s, payoff_typeValue, payoff_params):
    if payoff_typeValue == FXRainbowOptionTypes.CALL_ON_MINIMUM.value:
        k = payoff_params[0]
        payoff = np.maximum(np.min(s, axis=1) - k, 0.0)
    elif payoff_typeValue == FXRainbowOptionTypes.CALL_ON_MAXIMUM.value:
        k = payoff_params[0]
        payoff = np.maximum(np.max(s, axis=1) - k, 0.0)
    elif payoff_typeValue == FXRainbowOptionTypes.PUT_ON_MINIMUM.value:
        k = payoff_params[0]
        payoff = np.maximum(k - np.min(s, axis=1), 0.0)
    elif payoff_typeValue == FXRainbowOptionTypes.PUT_ON_MAXIMUM.value:
        k = payoff_params[0]
        payoff = np.maximum(k - np.max(s, axis=1), 0.0)
    elif payoff_typeValue == FXRainbowOptionTypes.CALL_ON_NTH.value:
        n = payoff_params[0]
        k = payoff_params[1]
        ssorted = np.sort(s)
        assetn = ssorted[:, -n]
        payoff = np.maximum(assetn - k, 0.0)
    elif payoff_typeValue == FXRainbowOptionTypes.PUT_ON_NTH.value:
        n = payoff_params[0]
        k = payoff_params[1]
        ssorted = np.sort(s)
        assetn = ssorted[:, -n]
        payoff = np.maximum(k - assetn, 0.0)
    else:
        raise FinError("Unknown payoff type")

    return payoff


###############################################################################


def value_mc_fast(
    t,
    stock_prices,
    discount_curve,
    dividend_yields,
    volatilities,
    betas,
    num_assets,
    payoff_type,
    payoff_params,
    num_paths=10000,
    seed=4242,
):

    np.random.seed(seed)
    df = discount_curve.df(t)
    r = -np.log(df) / t
    mus = r - dividend_yields

    num_time_steps = 2
    s_all = get_assets_paths_times(
        num_assets,
        num_paths,
        num_time_steps,
        t,
        mus,
        stock_prices,
        volatilities,
        betas,
        seed,
    )

    payoff = payoff_value(s_all, payoff_type.value, payoff_params)
    payoff = np.mean(payoff)
    v = payoff * np.exp(-r * t)
    return v


###############################################################################


class FXRainbowOption(FXOption):

    def __init__(
        self,
        expiry_dt: Date,
        payoff_type: FXRainbowOptionTypes,
        payoff_params: List[float],
        num_assets: int,
    ):

        check_argument_types(self.__init__, locals())

        self.validate_payoff(payoff_type, payoff_params, num_assets)

        self.expiry_dt = expiry_dt
        self.payoff_type = payoff_type
        self.payoff_params = payoff_params
        self.num_assets = num_assets

    ###########################################################################

    def validate(self, stock_prices, dividend_yields, volatilities, betas):

        if len(stock_prices) != self.num_assets:
            raise FinError(
                "Stock prices must be a vector of length "
                + str(self.num_assets)
            )

        if len(dividend_yields) != self.num_assets:
            raise FinError(
                "Dividend yields must be a vector of length "
                + str(self.num_assets)
            )

        if len(volatilities) != self.num_assets:
            raise FinError(
                "Volatilities must be a vector of length "
                + str(self.num_assets)
            )

        if len(betas) != self.num_assets:
            raise FinError(
                "Betas must be a vector of length " + str(self.num_assets)
            )

    ###########################################################################

    def validate_payoff(self, payoff_type, payoff_params, num_assets):

        num_params = 0

        if payoff_type == FXRainbowOptionTypes.CALL_ON_MINIMUM:
            num_params = 1
        elif payoff_type == FXRainbowOptionTypes.CALL_ON_MAXIMUM:
            num_params = 1
        elif payoff_type == FXRainbowOptionTypes.PUT_ON_MINIMUM:
            num_params = 1
        elif payoff_type == FXRainbowOptionTypes.PUT_ON_MAXIMUM:
            num_params = 1
        elif payoff_type == FXRainbowOptionTypes.CALL_ON_NTH:
            num_params = 2
        elif payoff_type == FXRainbowOptionTypes.PUT_ON_NTH:
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
            payoff_type == FXRainbowOptionTypes.CALL_ON_NTH
            or payoff_type == FXRainbowOptionTypes.PUT_ON_NTH
        ):
            n = payoff_params[0]
            if n < 1 or n > num_assets:
                raise FinError("Nth parameter must be 1 to " + str(num_assets))

    ###########################################################################

    def value(
        self,
        value_dt,
        stock_prices,
        domestic_curve,
        foreign_curve,
        volatilities,
        betas,
    ):

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self.expiry_dt:
            raise FinError("Valuation date after expiry date.")

        if domestic_curve.value_dt != value_dt:
            raise FinError(
                "Domestic Curve valuation date not same as option value date"
            )

        if foreign_curve.value_dt != value_dt:
            raise FinError(
                "Foreign Curve valuation date not same as option value date"
            )

        if self.num_assets != 2:
            raise FinError("Analytical results for two assets only.")

        if value_dt > self.expiry_dt:
            raise FinError("Value date after expiry date.")

        self.validate(stock_prices, foreign_curve, volatilities, betas)

        # Use result by Stulz (1982) given by Haug Page 211
        t = (self.expiry_dt - value_dt) / g_days_in_year

        df = domestic_curve.df(t)
        r = -np.log(df) / t

        q1 = foreign_curve[0]
        q2 = foreign_curve[1]
        rho = betas[0] ** 2
        s1 = stock_prices[0]
        s2 = stock_prices[1]
        b1 = r - q1
        b2 = r - q2
        v1 = volatilities[0]
        v2 = volatilities[1]
        k = self.payoff_params[0]

        v = np.sqrt(v1 * v1 + v2 * v2 - 2 * rho * v1 * v2)
        d = (np.log(s1 / s2) + (b1 - b2 + v * v / 2) * t) / v / np.sqrt(t)
        y1 = (np.log(s1 / k) + (b1 + v1 * v1 / 2) * t) / v1 / np.sqrt(t)
        y2 = (np.log(s2 / k) + (b2 + v2 * v2 / 2) * t) / v2 / np.sqrt(t)
        rho1 = (v1 - rho * v2) / v
        rho2 = (v2 - rho * v1) / v
        dq1 = np.exp(-q1 * t)
        dq2 = np.exp(-q2 * t)
        df = np.exp(-r * t)
        sqrtt = np.sqrt(t)

        if self.payoff_type == FXRainbowOptionTypes.CALL_ON_MAXIMUM:
            v = (
                s1 * dq1 * M(y1, d, rho1)
                + s2 * dq2 * M(y2, -d + v * sqrtt, rho2)
                - k
                * df
                * (1.0 - M(-y1 + v1 * np.sqrt(t), -y2 + v2 * sqrtt, rho))
            )
        elif self.payoff_type == FXRainbowOptionTypes.CALL_ON_MINIMUM:
            v = (
                s1 * dq1 * M(y1, -d, -rho1)
                + s2 * dq2 * M(y2, d - v * np.sqrt(t), -rho2)
                - k * df * M(y1 - v1 * np.sqrt(t), y2 - v2 * np.sqrt(t), rho)
            )
        elif self.payoff_type == FXRainbowOptionTypes.PUT_ON_MAXIMUM:
            cmax1 = s2 * dq2 + s1 * dq1 * N(d) - s2 * dq2 * N(d - v * sqrtt)
            cmax2 = (
                s1 * dq1 * M(y1, d, rho1)
                + s2 * dq2 * M(y2, -d + v * sqrtt, rho2)
                - k * df * (1.0 - M(-y1 + v1 * sqrtt, -y2 + v2 * sqrtt, rho))
            )
            v = k * df - cmax1 + cmax2
        elif self.payoff_type == FXRainbowOptionTypes.PUT_ON_MINIMUM:
            cmin1 = s1 * dq1 - s1 * dq1 * N(d) + s2 * dq2 * N(d - v * sqrtt)
            cmin2 = (
                s1 * dq1 * M(y1, -d, -rho1)
                + s2 * dq2 * M(y2, d - v * sqrtt, -rho2)
                - k * df * M(y1 - v1 * sqrtt, y2 - v2 * sqrtt, rho)
            )
            v = k * df - cmin1 + cmin2
        else:
            raise FinError("Unsupported FX Rainbow option type")

        return v

    ###########################################################################

    def value_mc(
        self,
        value_dt,
        expiry_dt,
        stock_prices,
        discount_curve,
        dividend_yields,
        volatilities,
        betas,
        num_paths=10000,
        seed=4242,
    ):

        self.validate(stock_prices, dividend_yields, volatilities, betas)

        if value_dt > expiry_dt:
            raise FinError("Value date after expiry date.")

        t = (self.expiry_dt - value_dt) / g_days_in_year

        v = value_mc_fast(
            t,
            stock_prices,
            discount_curve,
            dividend_yields,
            volatilities,
            betas,
            self.num_assets,
            self.payoff_type,
            self.payoff_params,
            num_paths,
            seed,
        )

        return v


###############################################################################
