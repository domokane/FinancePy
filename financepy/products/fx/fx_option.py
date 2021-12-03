##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from ...models.black_scholes import BlackScholes
from ...utils.global_vars import gDaysInYear

##########################################################################

bump = 1e-4

##########################################################################


class FXOption:
    """ Class that is used to perform perturbation risk for FX options. """

###############################################################################

    def delta(self, valuation_date, stock_price, discount_curve,
              dividend_curve, model):
        """ Calculate the option delta (FX rate sensitivity) by adding on a
        small bump and calculating the change in the option price. """

        v = self.value(valuation_date, stock_price, discount_curve, dividend_curve,
                       model)

        vBumped = self.value(valuation_date, stock_price + bump, discount_curve,
                             dividend_curve, model)

        if type(vBumped) is dict:
            delta = (vBumped['value'] - v['value']) / bump
        else:
            delta = (vBumped - v) / bump

        return delta

###############################################################################

    def gamma(self, valuation_date, stock_price, discount_curve, dividend_curve,
              model):
        """ Calculate the option gamma (delta sensitivity) by adding on a
        small bump and calculating the change in the option delta. """

        v = self.delta(valuation_date, stock_price, discount_curve, dividend_curve,
                       model)

        vBumpedDn = self.delta(valuation_date, stock_price + bump, discount_curve,
                               dividend_curve, model)

        vBumpedUp = self.delta(valuation_date, stock_price + bump, discount_curve,
                               dividend_curve, model)

        if type(v) is dict:
            num = (vBumpedUp['value'] - 2.0 * v['value'] + vBumpedDn['value'])
            gamma = num / bump / 2.0
        else:
            gamma = (vBumpedUp - 2.0 * v + vBumpedDn) / bump / 2.0

        return gamma

###############################################################################

    def vega(self, valuation_date, stock_price, discount_curve, dividend_curve, model):
        """ Calculate the option vega (volatility sensitivity) by adding on a
        small bump and calculating the change in the option price. """

        v = self.value(valuation_date, stock_price, discount_curve, dividend_curve,
                       model)

        vp = self.value(valuation_date, stock_price, discount_curve, dividend_curve,
                        BlackScholes(model._volatility + bump))

        if type(v) is dict:
            vega = (vp['value'] - v['value']) / bump
        else:
            vega = (vp - v) / bump

        return vega

###############################################################################

    def theta(self, valuation_date, stock_price, discount_curve, dividend_curve, model):
        """ Calculate the option theta (calendar time sensitivity) by moving
        forward one day and calculating the change in the option price. """

        v = self.value(valuation_date, stock_price, discount_curve,
                       dividend_curve, model)

        next_date = valuation_date.add_days(1)

        discount_curve._valuation_date = next_date
        dividend_curve._valuation_date = next_date

        vBumped = self.value(next_date, stock_price, discount_curve,
                             dividend_curve, model)

        bump = 1.0 / gDaysInYear

        if type(v) is dict:
            theta = (vBumped['value'] - v['value']) / bump
        else:
            theta = (vBumped - v) / bump

        discount_curve._valuation_date = valuation_date
        dividend_curve._valuation_date = valuation_date

        return theta

##############################################################################

    def rho(self, valuation_date, stock_price, discount_curve, dividend_curve, model):
        """ Calculate the option rho (interest rate sensitivity) by perturbing
        the discount curve and revaluing. """

        v = self.value(
            valuation_date,
            stock_price,
            discount_curve,
            dividend_curve,
            model)
        vBumped = self.value(
            valuation_date,
            stock_price,
            discount_curve.bump(bump),
            dividend_curve,
            model)

        if type(v) is dict:
            rho = (vBumped['value'] - v['value']) / bump
        else:
            rho = (vBumped - v) / bump

        return rho

##########################################################################
##########################################################################
