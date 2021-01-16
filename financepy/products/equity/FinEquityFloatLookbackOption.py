##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np


from ...finutils.FinMath import N
from ...finutils.FinGlobalVariables import gDaysInYear, gSmall
from ...finutils.FinError import FinError
from ...finutils.FinDate import FinDate

from ...models.FinGBMProcess import FinGBMProcess
from ...products.equity.FinEquityOption import FinEquityOption
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...market.curves.FinDiscountCurve import FinDiscountCurve
from ...finutils.FinGlobalTypes import FinOptionTypes

##########################################################################
# TODO: Attempt control variate adjustment to monte carlo
# TODO: Sobol for Monte Carlo
# TODO: TIGHTEN UP LIMIT FOR W FROM 100
# TODO: Vectorise the analytical pricing formula
##########################################################################


##########################################################################
# FLOAT STRIKE LOOKBACK CALL PAYS MAX(S(T)-SMIN,0)
# FLOAT STRIKE LOOKBACK PUT PAYS MAX(SMAX-S(T),0)
##########################################################################


class FinEquityFloatLookbackOption(FinEquityOption):
    ''' This is an equity option in which the strike of the option is not fixed
    but is set at expiry to equal the minimum stock price in the case of a call
    or the maximum stock price in the case of a put. In other words the buyer
    of the call gets to buy the asset at the lowest price over the period
    before expiry while the buyer of the put gets to sell the asset at the
    highest price before expiry. '''

    def __init__(self,
                 expiryDate: FinDate,
                 optionType: FinOptionTypes):
        ''' Create the FloatLookbackOption by specifying the expiry date and
        the option type. The strike is determined internally as the maximum or
        minimum of the stock price depending on whether it is a put or a call
        option. '''

        checkArgumentTypes(self.__init__, locals())

        if optionType != FinOptionTypes.EUROPEAN_CALL and optionType != FinOptionTypes.EUROPEAN_PUT:
            raise FinError("Option type must be EUROPEAN_CALL or EUROPEAN_PUT")

        self._expiryDate = expiryDate
        self._optionType = optionType

###############################################################################

    def value(self,
              valueDate: FinDate,
              stockPrice: float,
              discountCurve: FinDiscountCurve,
              dividendCurve: FinDiscountCurve,
              volatility: float,
              stockMinMax: float):
        ''' Valuation of the Floating Lookback option using Black-Scholes using
        the formulae derived by Goldman, Sosin and Gatto (1979). '''

        t = (self._expiryDate - valueDate) / gDaysInYear
        df = discountCurve.df(self._expiryDate)


        r = discountCurve.ccRate(self._expiryDate)
        q = dividendCurve.ccRate(self._expiryDate)

        v = volatility
        s0 = stockPrice
        smin = 0.0
        smax = 0.0

        if self._optionType == FinOptionTypes.EUROPEAN_CALL:
            smin = stockMinMax
            if smin > s0:
                raise FinError(
                    "Smin must be less than or equal to the stock price.")
        elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
            smax = stockMinMax
            if smax < s0:
                raise FinError(
                    "Smax must be greater than or equal to the stock price.")

        if abs(r - q) < gSmall:
            q = r + gSmall

        dq = np.exp(-q * t)
        df = np.exp(-r * t)
        b = r - q
        u = v * v / 2.0 / b
        w = 2.0 * b / v / v
        expbt = np.exp(b * t)

        # Taken from Haug Page 142
        if self._optionType == FinOptionTypes.EUROPEAN_CALL:

            a1 = (np.log(s0 / smin) + (b + (v**2) / 2.0) * t) / v / np.sqrt(t)
            a2 = a1 - v * np.sqrt(t)

            if smin == s0:
                term = N(-a1 + 2.0 * b * np.sqrt(t) / v) - expbt * N(-a1)
            elif s0 < smin and w < -100:
                term = - expbt * N(-a1)
            else:
                term = ((s0 / smin)**(-w))\
                    * N(-a1 + 2.0 * b * np.sqrt(t) / v) - expbt * N(-a1)

            v = s0 * dq * N(a1) - smin * df * N(a2) + s0 * df * u * term

        elif self._optionType == FinOptionTypes.EUROPEAN_PUT:

            b1 = (np.log(s0 / smax) + (b + (v**2) / 2.0) * t) / v / np.sqrt(t)
            b2 = b1 - v * np.sqrt(t)

            if smax == s0:
                term = -N(b1 - 2.0 * b * np.sqrt(t) / v) + expbt * N(b1)
            elif s0 < smax and w > 100:
                term = expbt * N(b1)
            else:
                term = (-(s0 / smax)**(-w)) * \
                    N(b1 - 2.0 * b * np.sqrt(t) / v) + expbt * N(b1)

            v = smax * df * N(-b2) - s0 * dq * N(-b1) + s0 * df * u * term

        else:
            raise FinError("Unknown lookback option type:" +
                           str(self._optionType))

        return v

###############################################################################

    def valueMC(self,
                valueDate: FinDate,
                stockPrice: float,
                discountCurve: FinDiscountCurve,
                dividendCurve: FinDiscountCurve,
                volatility: float,
                stockMinMax: float,
                numPaths: int = 10000,
                numStepsPerYear: int = 252,
                seed: int = 4242):
        ''' Monte Carlo valuation of a floating strike lookback option using a
        Black-Scholes model that assumes the stock follows a GBM process. '''

        t = (self._expiryDate - valueDate) / gDaysInYear
        numTimeSteps = int(t * numStepsPerYear)

        df = discountCurve.df(self._expiryDate)
        r = discountCurve.ccRate(self._expiryDate)
        q = dividendCurve.ccRate(self._expiryDate)
        mu = r - q

        optionType = self._optionType
        smin = 0.0
        smax = 0.0

        if self._optionType == FinOptionTypes.EUROPEAN_CALL:
            smin = stockMinMax
            if smin > stockPrice:
                raise FinError(
                    "Smin must be less than or equal to the stock price.")
        elif self._optionType == FinOptionTypes.EUROPEAN_PUT:
            smax = stockMinMax
            if smax < stockPrice:
                raise FinError(
                    "Smax must be greater than or equal to the stock price.")

        model = FinGBMProcess()
        Sall = model.getPaths(numPaths, numTimeSteps, t, mu, stockPrice,
                              volatility, seed)

        # Due to antithetics we have doubled the number of paths
        numPaths = 2 * numPaths
        payoff = np.zeros(numPaths)

        if optionType == FinOptionTypes.EUROPEAN_CALL:
            SMin = np.min(Sall, axis=1)
            SMin = np.minimum(SMin, smin)
            payoff = np.maximum(Sall[:, -1] - SMin, 0.0)
        elif optionType == FinOptionTypes.EUROPEAN_PUT:
            SMax = np.max(Sall, axis=1)
            SMax = np.maximum(SMax, smax)
            payoff = np.maximum(SMax - Sall[:, -1], 0.0)
        else:
            raise FinError("Unknown lookback option type:" + str(optionType))

        v = payoff.mean() * df
        return v

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("EXPIRY DATE", self._expiryDate)
        s += labelToString("OPTION TYPE", self._optionType, "")
        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
