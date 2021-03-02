##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from enum import Enum

from ...utils.global_variables import gDaysInYear
from ...models.black_scholes import FinModelBlackScholes
from ...market.curves.discount_curve import DiscountCurve
from ...utils.date import Date

###############################################################################

bump = 1e-4

###############################################################################

class FinEquityOptionModelTypes(Enum):
    BLACKSCHOLES = 1
    ANOTHER = 2

###############################################################################


class FinEquityOption(object):
    """ This class is a parent class for all option classes that require any
    perturbatory risk. """

###############################################################################

    def value(self,
              valuation_date: Date,
              stock_price: float,
              discount_curve: DiscountCurve,
              dividendYield: float,
              model):

        print("You should not be here!")
        return 0.0

###############################################################################

    def delta(self,
              valuation_date: Date,
              stock_price: float,
              discount_curve: DiscountCurve,
              dividendCurve: DiscountCurve,
              model):
        """ Calculation of option delta by perturbation of stock price and
        revaluation. """
        v = self.value(valuation_date, stock_price, discount_curve,
                       dividendCurve, model)

        vBumped = self.value(valuation_date, stock_price + bump, discount_curve,
                             dividendCurve, model)

        delta = (vBumped - v) / bump
        return delta

###############################################################################

    def gamma(self,
              valuation_date: Date,
              stock_price: float,
              discount_curve: DiscountCurve,
              dividendCurve: DiscountCurve,
              model):
        """ Calculation of option gamma by perturbation of stock price and
        revaluation. """

        v = self.value(valuation_date, stock_price, discount_curve,
                       dividendCurve, model)

        vBumpedDn = self.value(valuation_date, stock_price - bump, discount_curve,
                               dividendCurve, model)

        vBumpedUp = self.value(valuation_date, stock_price + bump, discount_curve,
                               dividendCurve, model)

        gamma = (vBumpedUp - 2.0 * v + vBumpedDn) / bump / bump
        return gamma

###############################################################################

    def vega(self,
             valuation_date: Date,
             stock_price: float,
             discount_curve: DiscountCurve,
             dividendCurve: DiscountCurve,
             model):
        """ Calculation of option vega by perturbing vol and revaluation. """

        v = self.value(valuation_date, stock_price, discount_curve,
                       dividendCurve, model)

        model = FinModelBlackScholes(model._volatility + bump)

        vBumped = self.value(valuation_date, stock_price, discount_curve,
                             dividendCurve, model)

        vega = (vBumped - v) / bump
        return vega

###############################################################################

    def theta(self,
              valuation_date: Date,
              stock_price: float,
              discount_curve: DiscountCurve,
              dividendCurve: DiscountCurve,
              model):
        """ Calculation of option theta by perturbing value date by one 
        calendar date (not a business date) and then doing revaluation and 
        calculating the difference divided by dt = 1 / gDaysInYear. """

        v = self.value(valuation_date, stock_price,
                       discount_curve,
                       dividendCurve, model)

        nextDate = valuation_date.addDays(1)

        # Need to do this carefully.

        discount_curve._valuation_date = nextDate
        bump = (nextDate - valuation_date) / gDaysInYear

        vBumped = self.value(nextDate, stock_price,
                             discount_curve,
                             dividendCurve, model)

        discount_curve._valuation_date = valuation_date
        theta = (vBumped - v) / bump
        return theta

###############################################################################

    def rho(self,
            valuation_date: Date,
            stock_price: float,
            discount_curve: DiscountCurve,
            dividendCurve: DiscountCurve,
            model):
        """ Calculation of option rho by perturbing interest rate and
        revaluation. """

        v = self.value(valuation_date, stock_price, discount_curve,
                       dividendCurve, model)

        vBumped = self.value(valuation_date, stock_price, discount_curve.bump(bump),
                             dividendCurve, model)

        rho = (vBumped - v) / bump
        return rho

###############################################################################
