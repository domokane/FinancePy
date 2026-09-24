##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from typing import Union

import numpy as np


from ...utils.error import FinError
from ...utils.global_types import OptionTypes
from ...utils.global_types import DigitalOptionTypes
from ...products.equity.equity_option import EquityOption
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.date import Date
from ...market.curves.discount_curve import DiscountCurve
from ...models.bs_digital_option import bs_digital_option_value
from ...utils.check_values import check_curve_dt
from ...utils.check_values import check_t_exp
from ...utils.helpers import option_years
from ...models.model import Model

########################################################################################


class EquityDigitalOption(EquityOption):
    """A EquityDigitalOption is an option in which the buyer receives some
    payment if the stock price has crossed a barrier ONLY at expiry and zero
    otherwise. There are two types: cash-or-nothing and the asset-or-nothing
    option. We do not care whether the stock price has crossed the barrier
    today, we only care about the barrier at option expiry. For a continuously-
    monitored barrier, use the EquityOneTouchOption class."""

    def __init__(
        self,
        expiry_dt: Date,
        barrier: float,
        call_put_type: OptionTypes,
        digital_type: DigitalOptionTypes,
    ) -> None:
        """Create the digital option by specifying the expiry date, the
        barrier price and the type of option which is either a EUROPEAN_CALL
        or a EUROPEAN_PUT or an AMERICAN_CALL or AMERICAN_PUT. There are two
        types of underlying - cash or nothing and asset or nothing."""

        check_argument_types(self.__init__, locals())

        if call_put_type not in [
            OptionTypes.EUROPEAN_CALL,
            OptionTypes.EUROPEAN_PUT,
        ]:
            raise FinError("Option type must be EUROPEAN CALL or PUT")

        self.expiry_dt = expiry_dt
        self.barrier = float(barrier)
        self.call_put_type = call_put_type
        self.digital_type = digital_type

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: Union[float, np.ndarray],
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
    ):
        """Digital Option valuation using the Black-Scholes model assuming a
        barrier at expiry. Handles both cash-or-nothing and asset-or-nothing
        options."""

        t_exp = option_years(value_dt, self.expiry_dt)
        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)

        v = bs_digital_option_value(
            stock_price, t_exp, self.barrier, r, q, model.volatility, self.call_put_type.value, self.digital_type.value
        )

        return v

    ###########################################################################

    def value_mc(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
        num_paths: int = 10000,
        seed: int = 4242,
    ):
        """Digital Option valuation using the Black-Scholes model and Monte
        Carlo simulation. Product assumes a barrier only at expiry. Monte Carlo
        handles both a cash-or-nothing and an asset-or-nothing option."""

        t_exp = check_t_exp(value_dt, self.expiry_dt)
        t_exp = max(t_exp, 1e-10)

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        np.random.seed(seed)

        df = discount_curve.df(self.expiry_dt)
        r = discount_curve.zero_rate_cc(self.expiry_dt)
        q = dividend_curve.zero_rate_cc(self.expiry_dt)
        mu = r - q

        volatility = model.volatility
        k = self.barrier
        sqrt_t_exp = np.sqrt(t_exp)

        # Use Antithetic variables
        g = np.random.normal(0.0, 1.0, size=num_paths)
        s = stock_price * np.exp((mu - volatility * volatility / 2.0) * t_exp)
        m = np.exp(g * sqrt_t_exp * volatility)

        s_1 = s * m
        s_2 = s / m

        payoff_a_1 = None
        payoff_a_2 = None

        if self.digital_type == DigitalOptionTypes.CASH_OR_NOTHING:
            if self.call_put_type == OptionTypes.EUROPEAN_CALL:
                payoff_a_1 = np.heaviside(s_1 - k, 0.0)
                payoff_a_2 = np.heaviside(s_2 - k, 0.0)
            elif self.call_put_type == OptionTypes.EUROPEAN_PUT:
                payoff_a_1 = np.heaviside(k - s_1, 0.0)
                payoff_a_2 = np.heaviside(k - s_2, 0.0)
        elif self.digital_type == DigitalOptionTypes.ASSET_OR_NOTHING:
            if self.call_put_type == OptionTypes.EUROPEAN_CALL:
                payoff_a_1 = s_1 * np.heaviside(s_1 - k, 0.0)
                payoff_a_2 = s_2 * np.heaviside(s_2 - k, 0.0)
            elif self.call_put_type == OptionTypes.EUROPEAN_PUT:
                payoff_a_1 = s_1 * np.heaviside(k - s_1, 0.0)
                payoff_a_2 = s_2 * np.heaviside(k - s_2, 0.0)

        payoff = np.mean(payoff_a_1) + np.mean(payoff_a_2)
        v = payoff * df / 2.0
        return v

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("BARRIER LEVEL", self.barrier)
        s += label_to_string("CALL-PUT TYPE", self.call_put_type)
        s += label_to_string("DIGITAL TYPE", self.digital_type, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
