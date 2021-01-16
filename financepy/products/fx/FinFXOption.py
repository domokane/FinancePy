##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from ...models.FinModelBlackScholes import FinModelBlackScholes
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinHelperFunctions import labelToString

##########################################################################

bump = 1e-4

##########################################################################


class FinFXOption(object):
    ''' Class that is used to perform perturbation risk for FX options. '''

###############################################################################

    def delta(self, valueDate, stockPrice, discountCurve,
              dividendCurve, model):
        ''' Calculate the option delta (FX rate sensitivity) by adding on a
        small bump and calculating the change in the option price. '''

        v = self.value(valueDate, stockPrice, discountCurve, dividendCurve,
                       model)

        vBumped = self.value(valueDate, stockPrice + bump, discountCurve,
                             dividendCurve, model)

        if type(vBumped) is dict:
            delta = (vBumped['value'] - v['value']) / bump
        else:
            delta = (vBumped - v) / bump

        return delta

###############################################################################

    def gamma(self, valueDate, stockPrice, discountCurve, dividendCurve,
              model):
        ''' Calculate the option gamma (delta sensitivity) by adding on a
        small bump and calculating the change in the option delta. '''

        v = self.delta(valueDate, stockPrice, discountCurve, dividendCurve,
            model)

        vBumpedDn = self.delta(valueDate, stockPrice + bump, discountCurve,
            dividendCurve, model)

        vBumpedUp = self.delta(valueDate, stockPrice + bump, discountCurve,
            dividendCurve, model)

        if type(v) is dict:
            num = (vBumpedUp['value'] - 2.0 * v['value'] + vBumpedDn['value'])
            gamma = num / bump / 2.0
        else:
            gamma = (vBumpedUp - 2.0 * v + vBumpedDn) / bump / 2.0

        return gamma

###############################################################################

    def vega(self, valueDate, stockPrice, discountCurve, dividendCurve, model):
        ''' Calculate the option vega (volatility sensitivity) by adding on a
        small bump and calculating the change in the option price. '''

        v = self.value(valueDate, stockPrice, discountCurve, dividendCurve,
                       model)

        vp = self.value(valueDate, stockPrice, discountCurve, dividendCurve,
                        FinModelBlackScholes(model._volatility + bump))

        if type(v) is dict:
            vega = (vp['value'] - v['value']) / bump
        else:
            vega = (vp - v) / bump

        return vega

###############################################################################

    def theta(self, valueDate, stockPrice, discountCurve, dividendCurve, model):
        ''' Calculate the option theta (calendar time sensitivity) by moving
        forward one day and calculating the change in the option price. '''

        v = self.value(valueDate, stockPrice, discountCurve,
                       dividendCurve, model)

        nextDate = valueDate.addDays(1)
        bump = 1.0 / gDaysInYear

        vBumped = self.value(nextDate, stockPrice, discountCurve,
                             dividendCurve, model)

        if type(v) is dict:
            theta = (vBumped['value'] - v['value']) / bump
        else:
            theta = (vBumped - v) / bump

        return theta

    def rho(self, valueDate, stockPrice, discountCurve, dividendCurve, model):
        ''' Calculate the option rho (interest rate sensitivity) by perturbing
        the discount curve and revaluing. '''

        v = self.value(
            valueDate,
            stockPrice,
            discountCurve,
            dividendCurve,
            model)
        vBumped = self.value(
            valueDate,
            stockPrice,
            discountCurve.bump(bump),
            dividendCurve,
            model)

        if type(v) is dict:
            rho = (vBumped['value'] - v['value']) / bump
        else:
            rho = (vBumped - v) / bump

        return rho

##########################################################################
##########################################################################
