##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from ...models.FinModelBlackScholes import FinModelBlackScholes
from ...utils.FinGlobalVariables import gDaysInYear
from ...utils.FinHelperFunctions import labelToString

##########################################################################

bump = 1e-4

##########################################################################


class FinFXOption(object):
    """ Class that is used to perform perturbation risk for FX options. """

###############################################################################

    def delta(self, valuation_date, stockPrice, discount_curve,
              dividendCurve, model):
        """ Calculate the option delta (FX rate sensitivity) by adding on a
        small bump and calculating the change in the option price. """

        v = self.value(valuation_date, stockPrice, discount_curve, dividendCurve,
                       model)

        vBumped = self.value(valuation_date, stockPrice + bump, discount_curve,
                             dividendCurve, model)

        if type(vBumped) is dict:
            delta = (vBumped['value'] - v['value']) / bump
        else:
            delta = (vBumped - v) / bump

        return delta

###############################################################################

    def gamma(self, valuation_date, stockPrice, discount_curve, dividendCurve,
              model):
        """ Calculate the option gamma (delta sensitivity) by adding on a
        small bump and calculating the change in the option delta. """

        v = self.delta(valuation_date, stockPrice, discount_curve, dividendCurve,
            model)

        vBumpedDn = self.delta(valuation_date, stockPrice + bump, discount_curve,
            dividendCurve, model)

        vBumpedUp = self.delta(valuation_date, stockPrice + bump, discount_curve,
            dividendCurve, model)

        if type(v) is dict:
            num = (vBumpedUp['value'] - 2.0 * v['value'] + vBumpedDn['value'])
            gamma = num / bump / 2.0
        else:
            gamma = (vBumpedUp - 2.0 * v + vBumpedDn) / bump / 2.0

        return gamma

###############################################################################

    def vega(self, valuation_date, stockPrice, discount_curve, dividendCurve, model):
        """ Calculate the option vega (volatility sensitivity) by adding on a
        small bump and calculating the change in the option price. """

        v = self.value(valuation_date, stockPrice, discount_curve, dividendCurve,
                       model)

        vp = self.value(valuation_date, stockPrice, discount_curve, dividendCurve,
                        FinModelBlackScholes(model._volatility + bump))

        if type(v) is dict:
            vega = (vp['value'] - v['value']) / bump
        else:
            vega = (vp - v) / bump

        return vega

###############################################################################

    def theta(self, valuation_date, stockPrice, discount_curve, dividendCurve, model):
        """ Calculate the option theta (calendar time sensitivity) by moving
        forward one day and calculating the change in the option price. """

        v = self.value(valuation_date, stockPrice, discount_curve,
                       dividendCurve, model)

        nextDate = valuation_date.addDays(1)
        bump = 1.0 / gDaysInYear

        vBumped = self.value(nextDate, stockPrice, discount_curve,
                             dividendCurve, model)

        if type(v) is dict:
            theta = (vBumped['value'] - v['value']) / bump
        else:
            theta = (vBumped - v) / bump

        return theta

    def rho(self, valuation_date, stockPrice, discount_curve, dividendCurve, model):
        """ Calculate the option rho (interest rate sensitivity) by perturbing
        the discount curve and revaluing. """

        v = self.value(
            valuation_date,
            stockPrice,
            discount_curve,
            dividendCurve,
            model)
        vBumped = self.value(
            valuation_date,
            stockPrice,
            discount_curve.bump(bump),
            dividendCurve,
            model)

        if type(v) is dict:
            rho = (vBumped['value'] - v['value']) / bump
        else:
            rho = (vBumped - v) / bump

        return rho

##########################################################################
##########################################################################
