##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO: Consider risk management
# TODO: Consider using Sobol
# TODO: Consider allowing weights on the individual basket assets
# TODO: Extend monte carlo to handle American options

import numpy as np

from ...finutils.FinGlobalVariables import gDaysInYear
from ...models.FinGBMProcess import FinGBMProcess

from ...finutils.FinError import FinError
from ...finutils.FinGlobalTypes import FinOptionTypes
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...finutils.FinHelperFunctions import _funcName
from ...finutils.FinDate import FinDate
from ...market.curves.FinDiscountCurve import FinDiscountCurve

from ...finutils.FinMath import N


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

    def _validate(self,
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

        if correlations.ndim != 2:
            raise FinError(
                "Correlation must be a 2D matrix ")

        if correlations.shape[0] != self._numAssets:
            raise FinError(
                "Correlation cols must have a length " + str(self._numAssets))

        if correlations.shape[1] != self._numAssets:
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
              valueDate: FinDate,
              stockPrices: np.ndarray,
              discountCurve: FinDiscountCurve,
              dividendCurves: (list),
              volatilities: np.ndarray,
              correlations: np.ndarray):
        ''' Basket valuation using a moment matching method to approximate the
        effective variance of the underlying basket value. This approach is
        able to handle a full rank correlation structure between the individual
        assets. '''

    # https://pdfs.semanticscholar.org/16ed/c0e804379e22ff36dcbab7e9bb06519faa43.pdf

        texp = (self._expiryDate - valueDate) / gDaysInYear

        if valueDate > self._expiryDate:
            raise FinError("Value date after expiry date.")

        qs = []
        for curve in dividendCurves:
            q = curve.ccRate(self._expiryDate)
            qs.append(q)

        v = volatilities
        s = stockPrices

        self._validate(stockPrices,
                       qs,
                       volatilities,
                       correlations)

        a = np.ones(self._numAssets) * (1.0 / self._numAssets)

        r = discountCurve.ccRate(self._expiryDate)        

        smean = 0.0
        for ia in range(0, self._numAssets):
            smean = smean + s[ia] * a[ia]

        lnS0k = np.log(smean / self._strikePrice)
        sqrtT = np.sqrt(texp)

        # Moment matching - starting with dividend
        qnum = 0.0
        qden = 0.0
        for ia in range(0, self._numAssets):
            qnum = qnum + a[ia] * s[ia] * np.exp(-qs[ia] * texp)
            qden = qden + a[ia] * s[ia]
        qhat = -np.log(qnum / qden) / texp

        # Moment matching - matching volatility
        vnum = 0.0
        for ia in range(0, self._numAssets):
            for ja in range(0, ia):
                rhoSigmaSigma = v[ia] * v[ja] * correlations[ia, ja]
                expTerm = (qs[ia] + qs[ja] - rhoSigmaSigma) * texp
                vnum = vnum + a[ia] * a[ja] * s[ia] * s[ja] * np.exp(-expTerm)

        vnum *= 2.0

        for ia in range(0, self._numAssets):
            rhoSigmaSigma = v[ia] ** 2
            expTerm = (2.0 * qs[ia] - rhoSigmaSigma) * texp
            vnum = vnum + ((a[ia] * s[ia]) ** 2) * np.exp(-expTerm)

        vhat2 = np.log(vnum / qnum / qnum) / texp

        den = np.sqrt(vhat2) * sqrtT
        mu = r - qhat
        d1 = (lnS0k + (mu + vhat2 / 2.0) * texp) / den
        d2 = (lnS0k + (mu - vhat2 / 2.0) * texp) / den

        if self._optionType == FinOptionTypes.EUROPEAN_CALL:
            v = smean * np.exp(-qhat * texp) * N(d1)
            v = v - self._strikePrice * np.exp(-r * texp) * N(d2)
        elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
            v = self._strikePrice * np.exp(-r * texp) * N(-d2)
            v = v - smean * np.exp(-qhat * texp) * N(-d1)
        else:
            raise FinError("Unknown option type")

        return v

###############################################################################

    def valueMC(self,
                valueDate: FinDate,
                stockPrices: np.ndarray,
                discountCurve: FinDiscountCurve,
                dividendCurves: (list),
                volatilities: np.ndarray,
                corrMatrix: np.ndarray,
                numPaths:int = 10000,
                seed:int = 4242):
        ''' Valuation of the EquityBasketOption using a Monte-Carlo simulation
        of stock prices assuming a GBM distribution. Cholesky decomposition is
        used to handle a full rank correlation structure between the individual
        assets. The numPaths and seed are pre-set to default values but can be
        overwritten. '''

        checkArgumentTypes(getattr(self, _funcName(), None), locals())

        if valueDate > self._expiryDate:
            raise FinError("Value date after expiry date.")

        texp = (self._expiryDate - valueDate) / gDaysInYear

        dividendYields = []
        for curve in dividendCurves:
            dq = curve.df(self._expiryDate)
            q = -np.log(dq) / texp
            dividendYields.append(q)

        self._validate(stockPrices,
                       dividendYields,
                       volatilities,
                       corrMatrix)

        numAssets = len(stockPrices)

        df = discountCurve.df(self._expiryDate)
        r = -np.log(df)/texp

        mus = r - dividendYields
        k = self._strikePrice

        numTimeSteps = 2

        model = FinGBMProcess()
        np.random.seed(seed)

        Sall = model.getPathsAssets(numAssets,
                                    numPaths,
                                    numTimeSteps,
                                    texp,
                                    mus,
                                    stockPrices,
                                    volatilities,
                                    corrMatrix,
                                    seed)

        if self._optionType == FinOptionTypes.EUROPEAN_CALL:
            payoff = np.maximum(np.mean(Sall, axis=1) - k, 0.0)
        elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
            payoff = np.maximum(k - np.mean(Sall, axis=1), 0.0)
        else:
            raise FinError("Unknown option type.")

        payoff = np.mean(payoff)
        v = payoff * np.exp(-r * texp)
        return v

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("EXPIRY DATE", self._expiryDate)
        s += labelToString("STRIKE PRICE", self._strikePrice)
        s += labelToString("OPTION TYPE", self._optionType)
        s += labelToString("NUM ASSETS", self._numAssets, "")
        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
