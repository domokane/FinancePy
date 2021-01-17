##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np
from enum import Enum


from ...finutils.FinGlobalVariables import gDaysInYear, gSmall
from ...finutils.FinError import FinError
from ...finutils.FinGlobalTypes import FinOptionTypes
from ...products.equity.FinEquityOption import FinEquityOption
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...finutils.FinDate import FinDate
from ...market.curves.FinDiscountCurve import FinDiscountCurve

from ...finutils.FinMath import NVect


###############################################################################


class FinDigitalOptionTypes(Enum):
    CASH_OR_NOTHING = 1,
    ASSET_OR_NOTHING = 2

###############################################################################


class FinEquityDigitalOption(FinEquityOption):
    ''' A FinEquityDigitalOption is an option in which the buyer receives some
    payment if the stock price has crossed a barrier ONLY at expiry and zero
    otherwise. There are two types: cash-or-nothing and the asset-or-nothing
    option. We do not care whether the stock price has crossed the barrier
    today, we only care about the barrier at option expiry. For a continuously-
    monitored barrier, use the FinEquityOneTouchOption class. '''

    def __init__(self,
                 expiryDate: FinDate,
                 barrierPrice: float,
                 optionType: FinOptionTypes,
                 underlyingType: FinDigitalOptionTypes):
        ''' Create the digital option by specifying the expiry date, the
        barrier price and the type of option which is either a EUROPEAN_CALL
        or a EUROPEAN_PUT or an AMERICAN_CALL or AMERICAN_PUT. There are two
        types of underlying - cash or nothing and asset or nothing. '''

        checkArgumentTypes(self.__init__, locals())

        if optionType != FinOptionTypes.EUROPEAN_CALL and optionType != FinOptionTypes.EUROPEAN_PUT:
            raise FinError("Option type must be EUROPEAN CALL or EUROPEAN PUT")

        self._expiryDate = expiryDate
        self._barrierPrice = float(barrierPrice)
        self._optionType = optionType
        self._underlyingType = underlyingType

###############################################################################

    def value(self,
              valueDate: FinDate,
              stockPrice: (float, np.ndarray),
              discountCurve: FinDiscountCurve,
              dividendCurve: FinDiscountCurve,
              model):
        ''' Digital Option valuation using the Black-Scholes model assuming a
        barrier at expiry. Handles both cash-or-nothing and asset-or-nothing
        options.'''

        if valueDate > self._expiryDate:
            raise FinError("Value date after expiry date.")

        t = (self._expiryDate - valueDate) / gDaysInYear
        t = max(t, 1e-6)

        S0 = stockPrice
        X = self._barrierPrice
        lnS0k = np.log(S0 / X)

        sqrtT = np.sqrt(t)

        df = discountCurve.df(self._expiryDate)
        r = -np.log(df)/t

        dq = dividendCurve.df(self._expiryDate)
        q = -np.log(dq)/t

        volatility = model._volatility

        if abs(volatility) < gSmall:
            volatility = gSmall

        d1 = (lnS0k + (r - q + volatility*volatility / 2.0) * t)
        d1 = d1 / volatility / sqrtT
        d2 = d1 - volatility * sqrtT

        if self._underlyingType == FinDigitalOptionTypes.CASH_OR_NOTHING:
            if self._optionType == FinOptionTypes.EUROPEAN_CALL:
                v = np.exp(-r * t) * NVect(d2)
            elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
                v = np.exp(-r * t) * NVect(-d2)
        elif self._underlyingType == FinDigitalOptionTypes.ASSET_OR_NOTHING:
            if self._optionType == FinOptionTypes.EUROPEAN_CALL:
                v = S0 * np.exp(-q * t) * NVect(d1)
            elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
                v = S0 * np.exp(-q * t) * NVect(-d1)
        else:
            raise FinError("Unknown underlying type.")

        return v

###############################################################################

    def valueMC(self,
                valueDate: FinDate,
                stockPrice: float,
                discountCurve: FinDiscountCurve,
                dividendCurve: FinDiscountCurve,
                model,
                numPaths: int = 10000,
                seed: int = 4242):
        ''' Digital Option valuation using the Black-Scholes model and Monte
        Carlo simulation. Product assumes a barrier only at expiry. Monte Carlo
        handles both a cash-or-nothing and an asset-or-nothing option.'''

        np.random.seed(seed)
        t = (self._expiryDate - valueDate) / gDaysInYear
        df = discountCurve.df(self._expiryDate)
        r = -np.log(df)/t

        dq = dividendCurve.df(self._expiryDate)
        q = -np.log(dq)/t

        volatility = model._volatility
        K = self._barrierPrice
        sqrtdt = np.sqrt(t)

        # Use Antithetic variables
        g = np.random.normal(0.0, 1.0, size=(1, numPaths))
        s = stockPrice * np.exp((r - q - volatility * volatility / 2.0) * t)
        m = np.exp(g * sqrtdt * volatility)

        s_1 = s * m
        s_2 = s / m

        if self._underlyingType == FinDigitalOptionTypes.CASH_OR_NOTHING:
            if self._optionType == FinOptionTypes.EUROPEAN_CALL:
                payoff_a_1 = np.heaviside(s_1 - K, 0.0)
                payoff_a_2 = np.heaviside(s_2 - K, 0.0)
            elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
                payoff_a_1 = np.heaviside(K - s_1, 0.0)
                payoff_a_2 = np.heaviside(K - s_2, 0.0)
        elif self._underlyingType == FinDigitalOptionTypes.ASSET_OR_NOTHING:
            if self._optionType == FinOptionTypes.EUROPEAN_CALL:
                payoff_a_1 = s_1 * np.heaviside(s_1 - K, 0.0)
                payoff_a_2 = s_2 * np.heaviside(s_2 - K, 0.0)
            elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
                payoff_a_1 = s_1 * np.heaviside(K - s_1, 0.0)
                payoff_a_2 = s_2 * np.heaviside(K - s_2, 0.0)

        payoff = np.mean(payoff_a_1) + np.mean(payoff_a_2)
        v = payoff * df / 2.0
        return v

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("EXPIRY DATE", self._expiryDate)
        s += labelToString("BARRIER LEVEL", self._barrierPrice)
        s += labelToString("OPTION TYPE", self._optionType)
        s += labelToString("UNDERLYING TYPE", self._underlyingType, "")
        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
