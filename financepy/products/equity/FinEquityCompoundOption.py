##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from math import exp, log, sqrt

import numpy as np
from scipy import optimize
from numba import njit


from ...utils.fin_math import N, phi2
from ...utils.global_variables import gDaysInYear, gSmall
from ...utils.FinError import FinError
from ...utils.FinGlobalTypes import FinOptionTypes

from ...products.equity.FinEquityOption import FinEquityOption
from ...products.equity.FinEquityVanillaOption import FinEquityVanillaOption
from ...market.curves.discount_curve_flat import DiscountCurve
from ...market.curves.discount_curve_flat import DiscountCurveFlat
from ...utils.helper_functions import labelToString, check_argument_types
from ...utils.date import Date

###############################################################################
# TODO: Vectorise pricer
# TODO: Monte Carlo pricer
###############################################################################


def _f(s0, *args):

    self = args[0]
    valuation_date = args[1]
    discount_curve = args[2]
    dividendCurve = args[3]
    model = args[4]
    value = args[5]

    if s0 <= 0.0:
        raise FinError("Unable to solve for stock price that fits K1")

    objFn = self.value(valuation_date,
                       s0,
                       discount_curve,
                       dividendCurve,
                       model) - value

    return objFn

###############################################################################


@njit(fastmath=True, cache=True, nogil=True)
def _valueOnce(stock_price,
               r,
               q,
               volatility,
               t1, t2,
               optionType1, optionType2,
               k1, k2,
               num_steps):

    if num_steps < 3:
        num_steps = 3

    # TODO EUROPEAN call-put works but AMERICAN call-put needs to be tested

    # Need equally spaced time intervals for a recombining tree
    # Downside is that we may not measure periods exactly
    dt = t2 / num_steps
    num_steps1 = int(t1 / dt)
    num_steps2 = num_steps - num_steps1
    dt1 = dt
    dt2 = dt

    # print("T1:",t1,"T2:",t2,"dt:",dt,"N1*dt",num_steps1*dt,"N*dt",num_steps*dt)

    # the number of nodes on the tree
    numNodes = int(0.5 * (num_steps + 1) * (num_steps + 2))
    optionValues = np.zeros(numNodes)

    u1 = np.exp(volatility * np.sqrt(dt))
    d1 = 1.0 / u1
    u2 = np.exp(volatility * np.sqrt(dt))
    d2 = 1.0 / u2

    probs = np.zeros(num_steps)
    periodDiscountFactors = np.zeros(num_steps)

    # store time independent information for later use in tree
    for iTime in range(0, num_steps1):
        a1 = np.exp((r - q) * dt1)
        probs[iTime] = (a1 - d1) / (u1 - d1)
        periodDiscountFactors[iTime] = np.exp(-r * dt1)

    for iTime in range(num_steps1, num_steps):
        a2 = np.exp((r - q) * dt2)
        probs[iTime] = (a2 - d2) / (u2 - d2)
        periodDiscountFactors[iTime] = np.exp(-r * dt2)

    stockValues = np.zeros(numNodes)
    stockValues[0] = stock_price
    sLow = stock_price

    for iTime in range(1, num_steps1 + 1):
        sLow *= d1
        s = sLow
        for iNode in range(0, iTime + 1):
            index = int(0.5 * iTime * (iTime + 1))
            stockValues[index + iNode] = s
            s = s * (u1 * u1)

    for iTime in range(num_steps1 + 1, num_steps + 1):
        sLow *= d2
        s = sLow
        for iNode in range(0, iTime + 1):
            index = int(0.5 * iTime * (iTime + 1))
            stockValues[index + iNode] = s
            s = s * (u2 * u2)

    # work backwards by first setting values at expiry date t2
    index = int(0.5 * num_steps * (num_steps + 1))

    for iNode in range(0, iTime + 1):
        s = stockValues[index + iNode]
        if optionType2 == FinOptionTypes.EUROPEAN_CALL\
           or optionType2 == FinOptionTypes.AMERICAN_CALL:
            optionValues[index + iNode] = max(s - k2, 0.0)
        elif optionType2 == FinOptionTypes.EUROPEAN_PUT\
           or optionType2 == FinOptionTypes.AMERICAN_PUT:
            optionValues[index + iNode] = max(k2 - s, 0.0)

    # begin backward steps from expiry at t2 to first expiry at time t1
    for iTime in range(num_steps - 1, num_steps1, -1):
        index = int(0.5 * iTime * (iTime + 1))
        for iNode in range(0, iTime + 1):
            s = stockValues[index + iNode]
            nextIndex = int(0.5 * (iTime + 1) * (iTime + 2))
            nextNodeDn = nextIndex + iNode
            nextNodeUp = nextIndex + iNode + 1
            vUp = optionValues[nextNodeUp]
            vDn = optionValues[nextNodeDn]
            futureExpectedValue = probs[iTime] * vUp
            futureExpectedValue += (1.0 - probs[iTime]) * vDn
            holdValue = periodDiscountFactors[iTime] * futureExpectedValue

            exerciseValue = 0.0 # NUMBA NEEDS HELP TO DETERMINE THE TYPE

            if optionType1 == FinOptionTypes.AMERICAN_CALL:
                exerciseValue = max(s - k2, 0.0)
            elif optionType1 == FinOptionTypes.AMERICAN_PUT:
                exerciseValue = max(k2 - s, 0.0)

            optionValues[index + iNode] = max(exerciseValue, holdValue)

    # Now do payoff at the end of the first expiry period at t1
    iTime = num_steps1
    index = int(0.5 * iTime * (iTime + 1))
    for iNode in range(0, iTime + 1):
        s = stockValues[index + iNode]
        nextIndex = int(0.5 * (iTime + 1) * (iTime + 2))
        nextNodeDn = nextIndex + iNode
        nextNodeUp = nextIndex + iNode + 1
        vUp = optionValues[nextNodeUp]
        vDn = optionValues[nextNodeDn]
        futureExpectedValue = probs[iTime] * vUp
        futureExpectedValue += (1.0 - probs[iTime]) * vDn
        holdValue = periodDiscountFactors[iTime] * futureExpectedValue

        if optionType1 == FinOptionTypes.EUROPEAN_CALL\
           or optionType1 == FinOptionTypes.AMERICAN_CALL:
            optionValues[index + iNode] = max(holdValue - k1, 0.0)
        elif optionType1 == FinOptionTypes.EUROPEAN_PUT\
           or optionType1 == FinOptionTypes.AMERICAN_PUT:
            optionValues[index + iNode] = max(k1 - holdValue, 0.0)

    # begin backward steps from t1 expiry to value date
    for iTime in range(num_steps1 - 1, -1, -1):
        index = int(0.5 * iTime * (iTime + 1))
        for iNode in range(0, iTime + 1):
            s = stockValues[index + iNode]
            nextIndex = int(0.5 * (iTime + 1) * (iTime + 2))
            nextNodeDn = nextIndex + iNode
            nextNodeUp = nextIndex + iNode + 1
            vUp = optionValues[nextNodeUp]
            vDn = optionValues[nextNodeDn]
            futureExpectedValue = probs[iTime] * vUp
            futureExpectedValue += (1.0 - probs[iTime]) * vDn
            holdValue = periodDiscountFactors[iTime] * futureExpectedValue

            exerciseValue = 0.0 # NUMBA NEEDS HELP TO DETERMINE THE TYPE

            if optionType1 == FinOptionTypes.AMERICAN_CALL:
                exerciseValue = max(holdValue - k1, 0.0)
            elif optionType1 == FinOptionTypes.AMERICAN_PUT:
                exerciseValue = max(k1 - holdValue, 0.0)

            optionValues[index + iNode] = max(exerciseValue, holdValue)

    verbose = False
    if verbose:
        print("num_steps1:", num_steps1)
        print("num_steps2:", num_steps2)
        print("u1:", u1, "u2:", u2)
        print("dfs", periodDiscountFactors)
        print("probs:", probs)
        print("s:", stockValues)
        print("v4:", optionValues)

    # We calculate all of the important Greeks in one go
    price = optionValues[0]
    delta = (optionValues[2] - optionValues[1]) / \
        (stockValues[2] - stockValues[1])
    deltaUp = (optionValues[5] - optionValues[4]) / \
        (stockValues[5] - stockValues[4])
    deltaDn = (optionValues[4] - optionValues[3]) / \
        (stockValues[4] - stockValues[3])
    gamma = (deltaUp - deltaDn) / (stockValues[2] - stockValues[1])
    theta = (optionValues[4] - optionValues[0]) / (2.0 * dt1)
    results = np.array([price, delta, gamma, theta])
    return results

###############################################################################


class FinEquityCompoundOption(FinEquityOption):
    """ A FinEquityCompoundOption is a compound option which allows the holder
    to either buy or sell another underlying option on a first expiry date that
    itself expires on a second expiry date. Both strikes are set at trade
    initiation. """

    def __init__(self,
                 cExpiryDate: Date,  # Compound Option expiry date
                 cOptionType: FinOptionTypes,  # Compound option type
                 cStrikePrice: float,  # Compound option strike
                 uExpiryDate: Date,  # Underlying option expiry date
                 uOptionType: FinOptionTypes,  # Underlying option type
                 uStrikePrice: float):  # Underlying option strike price
        """ Create the FinEquityCompoundOption by passing in the first and
        second expiry dates as well as the corresponding strike prices and
        option types. """

        check_argument_types(self.__init__, locals())

        if cExpiryDate > uExpiryDate:
            raise FinError("Compound expiry date must precede underlying expiry date")

        if cOptionType != FinOptionTypes.EUROPEAN_CALL and \
            cOptionType != FinOptionTypes.AMERICAN_CALL and \
            cOptionType != FinOptionTypes.EUROPEAN_PUT and \
            cOptionType != FinOptionTypes.AMERICAN_PUT:
                raise FinError("Compound option must be European or American call or put.")

        if uOptionType != FinOptionTypes.EUROPEAN_CALL and \
            uOptionType != FinOptionTypes.AMERICAN_CALL and \
            uOptionType != FinOptionTypes.EUROPEAN_PUT and \
            uOptionType != FinOptionTypes.AMERICAN_PUT:
                raise FinError("Underlying Option must be European or American call or put.")

        self._cExpiryDate = cExpiryDate
        self._cStrikePrice = float(cStrikePrice)
        self._cOptionType = cOptionType

        self._uExpiryDate = uExpiryDate
        self._uStrikePrice = float(uStrikePrice)
        self._uOptionType = uOptionType

###############################################################################

    def value(self,
              valuation_date: Date,
              stock_price: float,
              discount_curve: DiscountCurve,
              dividendCurve: DiscountCurve,
              model,
              num_steps: int = 200):
        """ Value the compound option using an analytical approach if it is
        entirely European style. Otherwise use a Tree approach to handle the
        early exercise. Solution by Geske (1977), Hodges and Selby (1987) and
        Rubinstein (1991). See also Haug page 132. """

        # If the option has any American feature then use the tree
        if self._cOptionType == FinOptionTypes.AMERICAN_CALL or\
            self._uOptionType == FinOptionTypes.AMERICAN_CALL or\
            self._cOptionType == FinOptionTypes.AMERICAN_PUT or\
                self._uOptionType == FinOptionTypes.AMERICAN_PUT:

            v = self._valueTree(valuation_date,
                                stock_price,
                                discount_curve,
                                dividendCurve,
                                model,
                                num_steps)

            return v[0]

        tc = (self._cExpiryDate - valuation_date) / gDaysInYear
        tu = (self._uExpiryDate - valuation_date) / gDaysInYear

        s0 = stock_price

        df = discount_curve.df(self._uExpiryDate)
        ru = -np.log(df)/tu

        # CHECK INTEREST RATES AND IF THERE SHOULD BE TWO RU AND RC ?????
        tc = np.maximum(tc, gSmall)
        tu = np.maximum(tc, tu)

        dq = dividendCurve.df(self._uExpiryDate)
        q = -np.log(dq)/tu

        v = np.maximum(model._volatility, gSmall)

        kc = self._cStrikePrice
        ku = self._uStrikePrice

        sstar = self._impliedStockPrice(s0,
                                        self._cExpiryDate,
                                        self._uExpiryDate,
                                        kc,
                                        ku,
                                        self._uOptionType,
                                        ru, q, model)

        a1 = (log(s0 / sstar) + (ru - q + (v**2) / 2.0) * tc) / v / sqrt(tc)
        a2 = a1 - v * sqrt(tc)
        b1 = (log(s0 / ku) + (ru - q + (v**2) / 2.0) * tu) / v / sqrt(tu)
        b2 = b1 - v * sqrt(tu)

        dqu = exp(-q * tu)
        dfc = exp(-ru * tc)
        dfu = exp(-ru * tu)
        c = sqrt(tc / tu)

        # Taken from Hull Page 532 (6th edition)

        CALL = FinOptionTypes.EUROPEAN_CALL
        PUT = FinOptionTypes.EUROPEAN_PUT

        if self._cOptionType == CALL and self._uOptionType == CALL:
            v = s0 * dqu * phi2(a1, b1, c) - ku * dfu * \
                phi2(a2, b2, c) - dfc * kc * N(a2)
        elif self._cOptionType == PUT and self._uOptionType == CALL:
            v = ku * dfu * phi2(-a2, b2, -c) - s0 * dqu * \
                phi2(-a1, b1, -c) + dfc * kc * N(-a2)
        elif self._cOptionType == CALL and self._uOptionType == PUT:
            v = ku * dfu * phi2(-a2, -b2, c) - s0 * dqu * \
                phi2(-a1, -b1, c) - dfc * kc * N(-a2)
        elif self._cOptionType == PUT and self._uOptionType == PUT:
            v = s0 * dqu * phi2(a1, -b1, -c) - ku * dfu * \
                phi2(a2, -b2, -c) + dfc * kc * N(a2)
        else:
            raise FinError("Unknown option type")

        return v

###############################################################################

    def _valueTree(self,
                   valuation_date,
                   stock_price,
                   discount_curve,
                   dividendCurve,
                   model,
                   num_steps=200):
        """ This function is called if the option has American features. """

        if valuation_date > self._cExpiryDate:
            raise FinError("Value date is after expiry date.")

        tc = (self._cExpiryDate - valuation_date) / gDaysInYear
        tu = (self._uExpiryDate - valuation_date) / gDaysInYear

        df = discount_curve.df(self._uExpiryDate)
        r = -np.log(df)/tu

        dq = dividendCurve.df(self._uExpiryDate)
        q = -np.log(dq)/tu

        r = discount_curve.zeroRate(self._uExpiryDate)

        volatility = model._volatility

        v1 = _valueOnce(stock_price,
                        r,
                        q,
                        volatility,
                        tc, tu,
                        self._cOptionType,
                        self._uOptionType,
                        self._cStrikePrice,
                        self._uStrikePrice,
                        num_steps)

        return v1

###############################################################################

    def _impliedStockPrice(self,
                           stock_price,
                           expiry_date1,
                           expiry_date2,
                           strikePrice1,
                           strikePrice2,
                           optionType2,
                           interestRate,
                           dividendYield,
                           model):

        option = FinEquityVanillaOption(expiry_date2, strikePrice2, optionType2)

        discount_curve = DiscountCurveFlat(expiry_date1, interestRate)
        dividendCurve = DiscountCurveFlat(expiry_date1, dividendYield)

        argtuple = (option, expiry_date1, discount_curve, dividendCurve,
                    model, strikePrice1)

        sigma = optimize.newton(_f, x0=stock_price, args=argtuple, tol=1e-8,
                                maxiter=50, fprime2=None)
        return sigma

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("CPD EXPIRY DATE", self._cExpiryDate)
        s += labelToString("CPD STRIKE PRICE", self._cStrikePrice)
        s += labelToString("CPD OPTION TYPE", self._cOptionType)
        s += labelToString("UND EXPIRY DATE", self._uExpiryDate)
        s += labelToString("UND STRIKE PRICE", self._uStrikePrice)
        s += labelToString("UND OPTION TYPE", self._uOptionType)
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
