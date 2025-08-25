##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from enum import Enum

from ...utils.global_vars import G_DAYS_IN_YEARS
from ...models.black_scholes import BlackScholes
from ...market.curves.discount_curve import DiscountCurve
from ...utils.date import Date

########################################################################################

BUMP = 1e-4

########################################################################################


class EquityOptionModelTypes(Enum):
    BLACKSCHOLES = 1
    ANOTHER = 2


########################################################################################


class EquityOption:
    """This class is a parent class for all equitu option classes that
    require any perturbatory risk."""

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_yield: float,
        model,
    ):

        print("You should not be here!")
        return 0.0

    ###########################################################################

    def delta(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model,
    ):
        """Calculation of option delta by perturbation of stock price and
        revaluation."""
        v = self.value(value_dt, stock_price, discount_curve, dividend_curve, model)

        v_bumped = self.value(
            value_dt, stock_price + BUMP, discount_curve, dividend_curve, model
        )

        delta = (v_bumped - v) / BUMP
        return delta

    ###########################################################################

    def gamma(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model,
    ):
        """Calculation of option gamma by perturbation of stock price and
        revaluation."""

        v = self.value(value_dt, stock_price, discount_curve, dividend_curve, model)

        v_bumped_dn = self.value(
            value_dt, stock_price - BUMP, discount_curve, dividend_curve, model
        )

        v_bumped_up = self.value(
            value_dt, stock_price + BUMP, discount_curve, dividend_curve, model
        )

        gamma = (v_bumped_up - 2.0 * v + v_bumped_dn) / BUMP / BUMP
        return gamma

    ###########################################################################

    def vega(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model,
    ):
        """Calculation of option vega by perturbing vol and revaluation."""

        # The bump should be 1% (0.01) not 1bp (0.0001)
        bump = 0.01

        v = self.value(value_dt, stock_price, discount_curve, dividend_curve, model)

        model = BlackScholes(model.volatility + bump)

        v_bumped = self.value(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )

        vega = v_bumped - v
        return vega

    ###########################################################################

    def vanna(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model,
    ):
        """Calculation of option vanna by perturbing delta with respect to the
        stock price volatility."""

        delta = self.delta(value_dt, stock_price, discount_curve, dividend_curve, model)

        model = BlackScholes(model.volatility + BUMP)

        delta_bumped = self.delta(
            value_dt, stock_price, discount_curve, dividend_curve, model
        )

        vanna = (delta_bumped - delta) / BUMP
        return vanna

    ###########################################################################

    def theta(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model,
    ):
        """Calculation of option theta by perturbing value date by one
        calendar date (not a business date) and then doing revaluation and
        calculating the difference divided by dt = 1 / G_DAYS_IN_YEARS."""

        v = self.value(value_dt, stock_price, discount_curve, dividend_curve, model)

        next_dt = value_dt.add_days(1)

        # Need to do this carefully. This is a bit hacky.
        discount_curve.value_dt = next_dt
        dividend_curve.value_dt = next_dt
        time_bump = (next_dt - value_dt) / G_DAYS_IN_YEARS

        v_bumped = self.value(
            next_dt, stock_price, discount_curve, dividend_curve, model
        )

        # restore valuation dates
        discount_curve.value_dt = value_dt
        dividend_curve.value_dt = value_dt

        theta = (v_bumped - v) / time_bump
        return theta

    ###########################################################################

    def rho(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model,
    ):
        """Calculation of option rho by perturbing interest rate and
        revaluation."""

        v = self.value(value_dt, stock_price, discount_curve, dividend_curve, model)

        v_bumped = self.value(
            value_dt,
            stock_price,
            discount_curve.bump(BUMP),
            dividend_curve,
            model,
        )

        rho = (v_bumped - v) / BUMP
        return rho


########################################################################################
