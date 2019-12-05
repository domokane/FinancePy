# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from enum import Enum

from ...finutils.FinGlobalVariables import gDaysInYear
from .FinEquityModelTypes import FinEquityModelBlackScholes
##########################################################################

bump = 1e-4


class FinOptionTypes(Enum):
    EUROPEAN_CALL = 1
    EUROPEAN_PUT = 2
    AMERICAN_CALL = 3
    AMERICAN_PUT = 4
    DIGITAL_CALL = 5
    DIGITAL_PUT = 6
    ASIAN_CALL = 7
    ASIAN_PUT = 8
    COMPOUND_CALL = 9
    COMPOUND_PUT = 10


class FinOptionModelTypes(Enum):
    BLACKSCHOLES = 1
    ANOTHER = 2

##########################################################################
##########################################################################


class FinOption(object):

    def delta(
            self,
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model):
        v = self.value(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)
        vBumped = self.value(
            valueDate,
            stockPrice + bump,
            discountCurve,
            dividendYield,
            model)

        if type(vBumped) is dict:
            delta = (vBumped['value'] - v['value']) / bump
        else:
            delta = (vBumped - v) / bump

        return delta

    def gamma(
            self,
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model):
        v = self.delta(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)
        vBumpedDn = self.delta(
            valueDate,
            stockPrice + bump,
            discountCurve,
            dividendYield,
            model)
        vBumpedUp = self.delta(
            valueDate,
            stockPrice + bump,
            discountCurve,
            dividendYield,
            model)

        if type(v) is dict:
            gamma = (vBumpedUp['value'] - 2.0 * v['value'] + vBumpedDn['value']) / bump / 2.0
        else:
            gamma = (vBumpedUp - 2.0 * v + vBumpedDn) / bump / 2.0

        return gamma

    def vega(
            self,
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model):
        v = self.value(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)
        vBumped = self.value(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            FinEquityModelBlackScholes(model._volatility + bump))

        if type(v) is dict:
            vega = (vBumped['value'] - v['value']) / bump
        else:
            vega = (vBumped - v) / bump

        return vega

    def theta(
            self,
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model):
        v = self.value(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)
        nextDate = valueDate.addDays(1)
        bump = 1.0 / gDaysInYear
        vBumped = self.value(
            nextDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)

        if type(v) is dict:
            theta = (vBumped['value'] - v['value']) / bump
        else:
            theta = (vBumped - v) / bump

        return theta

    def rho(
            self,
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model):

        v = self.value(
            valueDate,
            stockPrice,
            discountCurve,
            dividendYield,
            model)
        vBumped = self.value(
            valueDate,
            stockPrice,
            discountCurve.bump(bump),
            dividendYield,
            model)

        if type(v) is dict:
            rho = (vBumped['value'] - v['value']) / bump
        else:
            rho = (vBumped - v) / bump

        return rho

##########################################################################
##########################################################################
