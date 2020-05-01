# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 14:23:13 2016

@author: Dominic O'Kane
"""

import numpy as np
from math import log, exp, sqrt

from ...finutils.FinDate import FinDate
from ...finutils.FinMath import ONE_MILLION
from ...finutils.FinGlobalVariables import gDaysInYear
from .FinEquityOption import FinEquityModelBlackScholes
from .FinEquityVanillaOption import FinEquityVanillaOption, FinEquityOptionTypes

###############################################################################


class FinEquityVarianceSwap(object):
    ''' Class for managing an equity variance swap contract. '''

    def __init__(self,
                 startDate,
                 maturityDateOrTenor,
                 strikeVariance,
                 notional=ONE_MILLION,
                 payStrikeFlag=True):
        ''' Create variance swap contract. '''

        if type(startDate) != FinDate:
            raise ValueError("Settlement date must be a FinDate.")

        if type(maturityDateOrTenor) == FinDate:
            maturityDate = maturityDateOrTenor
        else:
            maturityDate = startDate.addTenor(maturityDateOrTenor)

        if startDate >= maturityDate:
            raise ValueError("Start date after or same as maturity date")

        self._startDate = startDate
        self._maturityDate = maturityDate
        self._strikeVariance = strikeVariance
        self._notional = notional
        self._payStrike = payStrikeFlag

        # Replication portfolio is stored
        self._numPutOptions = 0
        self._numCallOptions = 0
        self._putWts = []
        self._putStrikes = []
        self._callWts = []
        self._callStrikes = []

###############################################################################

    def value(self, 
              valuationDate, 
              realisedVar, 
              fairStrikeVar, 
              liborCurve):
        ''' Calculate the value of the variance swap based on the realised
        volatility to the valuation date, the forward looking implied
        volatility to the maturity date using the libor discount curve. '''

        t1 = (valuationDate - self._startDate) / gDaysInYear
        t2 = (self._maturityDate - self._startDate) / gDaysInYear

        expectedVariance = t1 * realisedVar/t2
        expectedVariance += (t2-t1) * fairStrikeVar / t2

        payoff = expectedVariance - self._strikeVariance

        df = liborCurve.df(self._maturityDate)
        v = payoff * self._notional * df
        return v

###############################################################################

    def fairStrikeApprox(self,
                         valuationDate,
                         fwdStockPrice,
                         strikes,
                         volatilities):
        ''' This is an approximation of the fair strike variance by Demeterfi
        et al. (1999) which assumes that sigma(K) = sigma(F) - b(K-F)/F where
        F is the forward stock price and sigma(F) is the ATM forward vol. '''

        f = fwdStockPrice

        # TODO Linear interpolation - to be revisited
        atmVol = np.interp(f, strikes, volatilities)
        tmat = (self._maturityDate - valuationDate)/gDaysInYear

        ''' Calculate the slope of the volatility curve by taking the end
        points in the volatilities and strikes to calculate the gradient.'''

        dvol = volatilities[-1] - volatilities[0]
        dK = strikes[-1] - strikes[0]
        b = f * dvol / dK
        var = (atmVol**2) * sqrt(1.0+3.0*tmat*(b**2))
        return var

###############################################################################

    def fairStrike(self,
                   valuationDate,
                   stockPrice,
                   dividendYield,
                   volatilityCurve,
                   numCallOptions,
                   numPutOptions,
                   strikeSpacing,
                   discountCurve,
                   useForward=True):
        ''' Calculate the implied variance according to the volatility surface
        using a static replication methodology with a specially weighted
        portfolio of put and call options across a range of strikes using the
        approximate method set out by Demeterfi et al. 1999. '''

        self._numPutOptions = numPutOptions
        self._numCallOptions = numCallOptions

        callType = FinEquityOptionTypes.EUROPEAN_CALL
        putType = FinEquityOptionTypes.EUROPEAN_PUT

        tmat = (self._maturityDate - valuationDate)/gDaysInYear

        df = discountCurve.df(tmat)
        r = - log(df)/tmat

        s0 = stockPrice
        g = exp(r*tmat)
        fwd = stockPrice * g

        # This fixes the centre strike of the replication options
        if useForward is True:
            sstar = fwd
        else:
            sstar = stockPrice

        ''' Replication argument from Demeterfi, Derman, Kamal and Zhou from
        Goldman Sachs Research notes March 1999. See Appendix A. This aim is
        to use calls and puts to approximate the payoff of a log contract '''

        minStrike = sstar - (numPutOptions+1) * strikeSpacing

        self._putWts = []
        self._putStrikes = []
        self._callWts = []
        self._callStrikes = []

        # if the lower strike is < 0 we go to as low as the strike spacing
        if minStrike < strikeSpacing:
            k = sstar
            klist = [sstar]
            while k >= strikeSpacing:
                k -= strikeSpacing
                klist.append(k)
            putK = np.array(klist)
            self._numPutOptions = len(putK) - 1
        else:
            putK = np.linspace(sstar, minStrike, numPutOptions+2)

        self._putStrikes = putK

        maxStrike = sstar + (numCallOptions+1) * strikeSpacing
        callK = np.linspace(sstar, maxStrike, numCallOptions+2)

        self._callStrikes = callK

        optionTotal = 2.0*(r*tmat - (s0*g/sstar-1.0) - log(sstar/s0))/tmat

        self._callWts = np.zeros(numCallOptions)
        self._putWts = np.zeros(numPutOptions)

        def f(x): return (2.0/tmat)*((x-sstar)/sstar-log(x/sstar))

        sumWts = 0.0
        for n in range(0, self._numPutOptions):
            kp = putK[n+1]
            k = putK[n]
            self._putWts[n] = (f(kp)-f(k))/(k-kp) - sumWts
            sumWts += self._putWts[n]

        sumWts = 0.0
        for n in range(0, self._numCallOptions):
            kp = callK[n+1]
            k = callK[n]
            self._callWts[n] = (f(kp)-f(k))/(kp-k) - sumWts
            sumWts += self._callWts[n]

        piPut = 0.0
        for n in range(0, numPutOptions):
            k = putK[n]
            vol = volatilityCurve.volatility(k)
            opt = FinEquityVanillaOption(self._maturityDate, k, putType)
            model = FinEquityModelBlackScholes(vol)
            v = opt.value(valuationDate, s0, discountCurve,
                          dividendYield, model)
            piPut += v * self._putWts[n]

        piCall = 0.0
        for n in range(0, numCallOptions):
            k = callK[n]
            vol = volatilityCurve.volatility(k)
            opt = FinEquityVanillaOption(self._maturityDate, k, callType)
            model = FinEquityModelBlackScholes(vol)
            v = opt.value(valuationDate, s0, discountCurve,
                          dividendYield, model)
            piCall += v * self._callWts[n]

        pi = piCall + piPut
        optionTotal += g * pi
        var = optionTotal

        return var

###############################################################################

    def realisedVariance(self, closePrices, useLogs=True):

        ''' Calculate the realised variance according to market standard
        calculations which can either use log or percentage returns.'''

        numObservations = len(closePrices)

        for i in range(0, numObservations):
            if closePrices[i] <= 0.0:
                raise ValueError("Stock prices must be greater than zero")

        cumX2 = 0.0

        if useLogs is True:
            for i in range(1, numObservations):
                x = log(closePrices[i]/closePrices[i-1])
                cumX2 += x*x
        else:
            for i in range(1, numObservations):
                x = (closePrices[i]-closePrices[i-1])/closePrices[i-1]
                cumX2 += x*x

        var = cumX2 * 252.0 / numObservations
        return var


###############################################################################

    def print(self):

        if self._numPutOptions == 0 and self._numCallOptions == 0:
            return

        print("TYPE", "STRIKE", "WEIGHT")
        for n in range(0, self._numPutOptions):
            k = self._putStrikes[n]
            wt = self._putWts[n]*self._notional
            print("PUT %7.2f %10.3f" % (k, wt))

        for n in range(0, self._numCallOptions):
            k = self._callStrikes[n]
            wt = self._callWts[n]*self._notional
            print("CALL %7.2f %10.3f" % (k, wt))

###############################################################################
