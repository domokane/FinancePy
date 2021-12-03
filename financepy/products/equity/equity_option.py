##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from enum import Enum

from ...utils.global_vars import gDaysInYear
from ...models.black_scholes import BlackScholes
from ...market.curves.discount_curve import DiscountCurve
from ...utils.date import Date

###############################################################################

bump = 1e-4

###############################################################################


class EquityOptionModelTypes(Enum):
    BLACKSCHOLES = 1
    ANOTHER = 2

###############################################################################


class EquityOption:
    """ This class is a parent class for all option classes that require any
    perturbatory risk. """

###############################################################################

    def value(self,
              valuation_date: Date,
              stock_price: float,
              discount_curve: DiscountCurve,
              dividend_yield: float,
              model):

        print("You should not be here!")
        return 0.0

###############################################################################

    def delta(self,
              valuation_date: Date,
              stock_price: float,
              discount_curve: DiscountCurve,
              dividend_curve: DiscountCurve,
              model):
        """ Calculation of option delta by perturbation of stock price and
        revaluation. """
        v = self.value(valuation_date, stock_price, discount_curve,
                       dividend_curve, model)

        vBumped = self.value(valuation_date, stock_price + bump, discount_curve,
                             dividend_curve, model)

        delta = (vBumped - v) / bump
        return delta

###############################################################################

    def gamma(self,
              valuation_date: Date,
              stock_price: float,
              discount_curve: DiscountCurve,
              dividend_curve: DiscountCurve,
              model):
        """ Calculation of option gamma by perturbation of stock price and
        revaluation. """

        v = self.value(valuation_date, stock_price, discount_curve,
                       dividend_curve, model)

        vBumpedDn = self.value(valuation_date, stock_price - bump, discount_curve,
                               dividend_curve, model)

        vBumpedUp = self.value(valuation_date, stock_price + bump, discount_curve,
                               dividend_curve, model)

        gamma = (vBumpedUp - 2.0 * v + vBumpedDn) / bump / bump
        return gamma

###############################################################################

    def vega(self,
             valuation_date: Date,
             stock_price: float,
             discount_curve: DiscountCurve,
             dividend_curve: DiscountCurve,
             model):
        """ Calculation of option vega by perturbing vol and revaluation. """

        v = self.value(valuation_date, stock_price, discount_curve,
                       dividend_curve, model)

        model = BlackScholes(model._volatility + bump)

        vBumped = self.value(valuation_date, stock_price, discount_curve,
                             dividend_curve, model)

        vega = (vBumped - v) / bump
        return vega

##############################################################################

    def vanna(self,
              valuation_date: Date,
              stock_price: float,
              discount_curve: DiscountCurve,
              dividend_curve: DiscountCurve,
              model):
        """ Calculation of option vanna by perturbing delta with respect to the
        stock price volatility. """

        delta = self.delta(valuation_date,
                           stock_price,
                           discount_curve,
                           dividend_curve,
                           model)

        model = BlackScholes(model._volatility + bump)

        deltaBumped = self.delta(valuation_date,
                                 stock_price,
                                 discount_curve,
                                 dividend_curve,
                                 model)

        vanna = (deltaBumped - delta) / bump
        return vanna

###############################################################################

    def theta(self,
              valuation_date: Date,
              stock_price: float,
              discount_curve: DiscountCurve,
              dividend_curve: DiscountCurve,
              model):
        """ Calculation of option theta by perturbing value date by one 
        calendar date (not a business date) and then doing revaluation and 
        calculating the difference divided by dt = 1 / gDaysInYear. """

        v = self.value(valuation_date, stock_price,
                       discount_curve,
                       dividend_curve, model)

        next_date = valuation_date.add_days(1)

        # Need to do this carefully. This is a bit hacky.
        discount_curve._valuation_date = next_date
        dividend_curve._valuation_date = next_date
        bump = (next_date - valuation_date) / gDaysInYear

        vBumped = self.value(next_date, stock_price,
                             discount_curve,
                             dividend_curve, model)

        # restore valuation dates
        discount_curve._valuation_date = valuation_date
        dividend_curve._valuation_date = valuation_date

        theta = (vBumped - v) / bump
        return theta

###############################################################################

    def rho(self,
            valuation_date: Date,
            stock_price: float,
            discount_curve: DiscountCurve,
            dividend_curve: DiscountCurve,
            model):
        """ Calculation of option rho by perturbing interest rate and
        revaluation. """

        v = self.value(valuation_date, stock_price, discount_curve,
                       dividend_curve, model)

        vBumped = self.value(valuation_date, stock_price, discount_curve.bump(bump),
                             dividend_curve, model)

        rho = (vBumped - v) / bump
        return rho

###############################################################################
