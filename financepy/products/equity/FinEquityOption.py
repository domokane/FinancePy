##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from enum import Enum

from ...finutils.FinGlobalVariables import gDaysInYear
from .FinEquityModelTypes import FinEquityModelBlackScholes
from ...market.curves.FinDiscountCurve import FinDiscountCurve
from ...finutils.FinDate import FinDate

###############################################################################

bump = 1e-4


class FinEquityOptionModelTypes(Enum):
    BLACKSCHOLES = 1
    ANOTHER = 2

###############################################################################


class FinEquityOption(object):
    ''' This class is a parent class for all option classes that require any
    perturbatory risk. '''

###############################################################################

    def value(self,
              valueDate: FinDate,
              stockPrice: float,
              discountCurve: FinDiscountCurve,
              dividendYield: float,
              model):

        print("You should not be here!")
        return 0.0

###############################################################################

    def delta(self,
              valueDate: FinDate,
              stockPrice: float,
              discountCurve: FinDiscountCurve,
              dividendYield: float,
              model):
        ''' Calculation of option delta by perturbation of stock proce and
        revaluation. '''
        v = self.value(valueDate, stockPrice, discountCurve,
                       dividendYield, model)

        vBumped = self.value(valueDate, stockPrice + bump, discountCurve,
                             dividendYield, model)

        delta = (vBumped - v) / bump
        return delta

###############################################################################

    def gamma(self,
              valueDate: FinDate,
              stockPrice: float,
              discountCurve: FinDiscountCurve,
              dividendYield: float,
              model):
        ''' Calculation of option gamma by perturbation of stock price and
        revaluation. '''

        v = self.value(valueDate, stockPrice, discountCurve,
                       dividendYield, model)

        vBumpedDn = self.value(valueDate, stockPrice - bump, discountCurve,
                               dividendYield, model)

        vBumpedUp = self.value(valueDate, stockPrice + bump, discountCurve,
                               dividendYield, model)

        gamma = (vBumpedUp - 2.0 * v + vBumpedDn) / bump / bump
        return gamma

###############################################################################

    def vega(self,
             valueDate: FinDate,
             stockPrice: float,
             discountCurve: FinDiscountCurve,
             dividendYield: float,
             model):
        ''' Calculation of option vega by perturbing vol and revaluation. '''

        v = self.value(valueDate, stockPrice, discountCurve,
                       dividendYield, model)

        model = FinEquityModelBlackScholes(model._volatility + bump)

        vBumped = self.value(valueDate, stockPrice, discountCurve,
                             dividendYield, model)

        vega = (vBumped - v) / bump
        return vega

###############################################################################

    def theta(self,
              valueDate: FinDate,
              stockPrice: float,
              discountCurve: FinDiscountCurve,
              dividendYield: float,
              model):
        ''' Calculation of option theta by perturbing value date and
        revaluation. '''

        v = self.value(valueDate, stockPrice, discountCurve,
                       dividendYield, model)

        nextDate = valueDate.addDays(1)
        # Need to do this carefully.
        discountCurve._valuationDate = nextDate
        bump = (nextDate - valueDate) / gDaysInYear

        vBumped = self.value(nextDate, stockPrice, discountCurve,
                             dividendYield, model)

        discountCurve._valuationDate = valueDate
        theta = (vBumped - v) / bump
        return theta

###############################################################################

    def rho(self,
            valueDate: FinDate,
            stockPrice: float,
            discountCurve: FinDiscountCurve,
            dividendYield: float,
            model):
        ''' Calculation of option rho by perturbing interest rate and
        revaluation. '''

        v = self.value(valueDate, stockPrice, discountCurve,
                       dividendYield, model)

        vBumped = self.value(valueDate, stockPrice, discountCurve.bump(bump),
                             dividendYield, model)

        rho = (vBumped - v) / bump
        return rho

###############################################################################
