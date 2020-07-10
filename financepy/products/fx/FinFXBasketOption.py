##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from math import exp, log, sqrt
import numpy as np


from ...finutils.FinHelperFunctions import labelToString

from ...finutils.FinMath import N
from ...finutils.FinGlobalVariables import gDaysInYear
from ...models.FinGBMProcess import FinGBMProcess

###############################################################################

from ...products.fx.FinFXOption import FinFXOption
from ...finutils.FinOptionTypes import FinOptionTypes
from ...finutils.FinHelperFunctions import checkArgumentTypes
from ...finutils.FinDate import FinDate
from ...finutils.FinError import FinError


class FinFXBasketOption(FinFXOption):
    ''' Class to manage FX Basket Option which is an option on a portfolio of
    FX rates. '''

    def __init__(self,
                 expiryDate: FinDate,
                 strikePrice: float,
                 optionType: FinOptionTypes,
                 numAssets: int,
                 notional: float = 1.0):
        ''' Create FX Basket Option with expiry date, strike price, option
        type, number of assets and notional. '''

        checkArgumentTypes(self.__init__, locals())

        self._expiryDate = expiryDate
        self._strikePrice = float(strikePrice)
        self._optionType = optionType
        self._numAssets = numAssets
        self._notional = notional

###############################################################################

    def validate(self,
                 stockPrices,
                 dividendYields,
                 volatilities,
                 betas):
        ''' Check that there is an input for each asset in the basket. '''

        if len(stockPrices) != self._numAssets:
            raise FinError("Stock prices must be a vector of length "
                           + str(self._numAssets))

        if len(dividendYields) != self._numAssets:
            raise FinError("Dividend yields must be a vector of length "
                           + str(self._numAssets))

        if len(volatilities) != self._numAssets:
            raise FinError("Volatilities must be a vector of length "
                           + str(self._numAssets))

        if len(betas) != self._numAssets:
            raise FinError("Betas must be a vector of length "
                           + str(self._numAssets))

###############################################################################

    def value(self,
              valueDate,
              stockPrices,
              discountCurve,
              dividendYields,
              volatilities,
              betas):
        ''' Value an FX Basket Option using Black-Scholes closed-form model
        which takes into account mean and variance of underlying. '''

        if valueDate > self._expiryDate:
            raise FinError("Value date after expiry date.")

        self.validate(stockPrices,
                      dividendYields,
                      volatilities,
                      betas)

        q = dividendYields
        v = volatilities
        s = stockPrices

        a = np.ones(self._numAssets) * (1.0 / self._numAssets)

        t = (self._expiryDate - valueDate) / gDaysInYear

        df = discountCurve._df(t)
        r = -np.log(df)/t

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
                rhoSigmaSigma = v[ia] * v[ja] * betas[ia] * betas[ja]
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
                domDiscountCurve,
                forDiscountCurve,
                volatilities,
                betas,
                numPaths=10000,
                seed=4242):
        ''' Value the FX Basket Option using Monte Carlo. '''

        if valueDate > self._expiryDate:
            raise FinError("Value date after expiry date.")

        numAssets = len(stockPrices)
        t = (self._expiryDate - valueDate) / gDaysInYear
        df = domDiscountCurve._df(t)
        rd = -log(df)/t
        dq = forDiscountCurve._df(t)
        rf = -log(dq)/t
        mus = rd - rf
        K = self._strikePrice

        # We jump to the expiry date as this option is not path dependent
        numTimeSteps = 2

        model = FinGBMProcess()
        np.random.seed(seed)

        Sall = model.getPathsAssets(numAssets,
                                    numPaths,
                                    numTimeSteps,
                                    t,
                                    mus,
                                    stockPrices,
                                    volatilities,
                                    betas,
                                    seed)

        if self._optionType == FinOptionTypes.EUROPEAN_CALL:
            payoff = np.maximum(np.mean(Sall, axis=1) - K, 0)
        elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
            payoff = np.maximum(K - np.mean(Sall, axis=1), 0)
        else:
            raise FinError("Unknown option type.")

        payoff = np.mean(payoff)
        v = payoff * exp(-rd * t)
        return v

###############################################################################
