##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from math import exp, log, sqrt
import numpy as np


from ...finutils.FinMath import N
from ...finutils.FinGlobalVariables import gDaysInYear
from ...models.FinGBMProcess import FinGBMProcess

##########################################################################

from ...finutils.FinError import FinError
from ...finutils.FinOptionTypes import FinOptionTypes
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...finutils.FinDate import FinDate

###############################################################################


class FinEquityBasketOption():
    ''' A FinEquityBasketOption is a contract to buy a put or a call option on
    an equally weighted portfolio of different stocks, each with its own price,
    volatility and dividend yield. An analytical and monte-carlo pricing model
    have been implemented for a European style option. '''

    def __init__(self,
                 expiryDate: FinDate,
                 strikePrice: float,
                 optionType: FinOptionTypes,
                 numAssets: int):
        ''' Define the FinEquityBasket option by specifying its expiry date,
        its strike price, whether it is a put or call, and the number of
        underlying stocks in the basket. '''

        checkArgumentTypes(self.__init__, locals())

        self._expiryDate = expiryDate
        self._strikePrice = float(strikePrice)
        self._optionType = optionType
        self._numAssets = numAssets

###############################################################################

    def validate(self,
                 stockPrices,
                 dividendYields,
                 volatilities,
                 correlations):

        if len(stockPrices) != self._numAssets:
            raise FinError(
                "Stock prices must have a length " + str(self._numAssets))

        if len(dividendYields) != self._numAssets:
            raise FinError(
                "Dividend yields must have a length " + str(self._numAssets))

        if len(volatilities) != self._numAssets:
            raise FinError(
                "Volatilities must have a length " + str(self._numAssets))

        if len(correlations) != self._numAssets:
            raise FinError(
                "Correlation cols must have a length " + str(self._numAssets))

        if len(correlations[0]) != self._numAssets:
            raise FinError(
                "correlation rows must have a length " + str(self._numAssets))

        for i in range(0, self._numAssets):
            if correlations[i, i] != 1.0:
                raise FinError("Corr matrix must have 1.0 on the diagonal")

            for j in range(0, i):
                if abs(correlations[i, j]) > 1.0:
                    raise FinError("Correlations must be [-1, +1]")

                if abs(correlations[j, i]) > 1.0:
                    raise FinError("Correlations must be [-1, +1]")

                if correlations[i, j] != correlations[j, i]:
                    raise FinError("Correlation matrix must be symmetric")

###############################################################################

    def value(self,
              valueDate,
              stockPrices,
              discountCurve,
              dividendYields,
              volatilities,
              correlations):
        ''' Basket valuation using a moment matching method to approximate the
        effective variance of the underlying basket value. '''

    # https://pdfs.semanticscholar.org/16ed/c0e804379e22ff36dcbab7e9bb06519faa43.pdf

        if valueDate > self._expiryDate:
            raise FinError("Value date after expiry date.")

        self.validate(stockPrices,
                      dividendYields,
                      volatilities,
                      correlations)

        q = dividendYields
        v = volatilities
        s = stockPrices

        a = np.ones(self._numAssets) * (1.0 / self._numAssets)
        t = (self._expiryDate - valueDate) / gDaysInYear
        df = discountCurve.df(self._expiryDate)
        r = -np.log(df) / t

        smean = 0.0
        for ia in range(0, self._numAssets):
            smean = smean + s[ia] * a[ia]

        lnS0k = log(smean / self._strikePrice)
        sqrtT = sqrt(t)

        # Moment matching - starting with dividend
        qnum = 0.0
        qden = 0.0
        for ia in range(0, self._numAssets):
            qnum = qnum + a[ia] * s[ia] * exp(-q[ia] * t)
            qden = qden + a[ia] * s[ia]
        qhat = -log(qnum / qden) / t

        # Moment matching - matching volatility
        vnum = 0.0
        for ia in range(0, self._numAssets):
            for ja in range(0, ia):
                rhoSigmaSigma = v[ia] * v[ja] * correlations[ia, ja]
                expTerm = (q[ia] + q[ja] - rhoSigmaSigma) * t
                vnum = vnum + a[ia] * a[ja] * s[ia] * s[ja] * exp(-expTerm)

        vnum *= 2.0

        for ia in range(0, self._numAssets):
            rhoSigmaSigma = v[ia] ** 2
            expTerm = (2.0 * q[ia] - rhoSigmaSigma) * t
            vnum = vnum + ((a[ia] * s[ia]) ** 2) * exp(-expTerm)

        vhat2 = log(vnum / qnum / qnum) / t

        den = sqrt(vhat2) * sqrtT
        mu = r - qhat
        d1 = (lnS0k + (mu + vhat2 / 2.0) * t) / den
        d2 = (lnS0k + (mu - vhat2 / 2.0) * t) / den

        if self._optionType == FinOptionTypes.EUROPEAN_CALL:
            v = smean * exp(-qhat * t) * N(d1)
            v = v - self._strikePrice * exp(-r * t) * N(d2)
        elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
            v = self._strikePrice * exp(-r * t) * N(-d2)
            v = v - smean * exp(-qhat * t) * N(-d1)
        else:
            raise FinError("Unknown option type")

        return v

###############################################################################

    def valueMC(self,
                valueDate,
                stockPrices,
                discountCurve,
                dividendYields,
                volatilities,
                corrMatrix,
                numPaths=10000,
                seed=4242):

        if valueDate > self._expiryDate:
            raise FinError("Value date after expiry date.")

        self.validate(stockPrices,
                      dividendYields,
                      volatilities,
                      corrMatrix)

        numAssets = len(stockPrices)

        t = (self._expiryDate - valueDate) / gDaysInYear
        r = discountCurve.zeroRate(self._expiryDate)
        mus = r - dividendYields
        k = self._strikePrice

        numTimeSteps = 2

        model = FinGBMProcess()
        np.random.seed(seed)

        Sall = model.getPathsAssets(
            numAssets,
            numPaths,
            numTimeSteps,
            t,
            mus,
            stockPrices,
            volatilities,
            corrMatrix,
            seed)

        if self._optionType == FinOptionTypes.EUROPEAN_CALL:
            payoff = np.maximum(np.mean(Sall, axis=1) - k, 0)
        elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
            payoff = np.maximum(k - np.mean(Sall, axis=1), 0)
        else:
            raise FinError("Unknown option type.")

        payoff = np.mean(payoff)
        v = payoff * exp(-r * t)
        return v

###############################################################################

    def __repr__(self):
        s = labelToString("EXPIRY DATE", self._expiryDate)
        s += labelToString("STRIKE PRICE", self._strikePrice)
        s += labelToString("OPTION TYPE", self._optionType)
        s += labelToString("NUM ASSETS", self._numAssets, "")
        return s

###############################################################################

    def print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
