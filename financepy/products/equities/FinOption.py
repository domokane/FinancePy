# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from enum import Enum

from ...finutils.FinGlobalVariables import gDaysInYear

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
            interestRate,
            dividendYield,
            volatility):
        v = self.value(
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            volatility)
        vBumped = self.value(
            valueDate,
            stockPrice + bump,
            interestRate,
            dividendYield,
            volatility)
        delta = (vBumped - v) / bump
        return delta

    def gamma(
            self,
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            volatility):
        v = self.delta(
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            volatility)
        vBumped = self.delta(
            valueDate,
            stockPrice + bump,
            interestRate,
            dividendYield,
            volatility)
        gamma = (vBumped - v) / bump
        return gamma

    def vega(
            self,
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            volatility):
        v = self.value(
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            volatility)
        vBumped = self.value(
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            volatility + bump)
        vega = (vBumped - v) / bump
        return vega

    def theta(
            self,
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            volatility):
        v = self.value(
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            volatility)
        nextDate = valueDate.addDays(1)
        bump = 1.0 / gDaysInYear
        vBumped = self.value(
            nextDate,
            stockPrice,
            interestRate,
            dividendYield,
            volatility)
        theta = (vBumped - v) / bump
        return theta

    def rho(
            self,
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            volatility):
        v = self.value(
            valueDate,
            stockPrice,
            interestRate,
            dividendYield,
            volatility)
        vBumped = self.value(
            valueDate,
            stockPrice,
            interestRate + bump,
            dividendYield,
            volatility)
        rho = (vBumped - v) / bump
        return rho

##########################################################################
##########################################################################
