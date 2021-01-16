##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np
from scipy import optimize

from ...finutils.FinMath import M
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinGlobalVariables import gSmall
from ...finutils.FinError import FinError

from ...products.equity.FinEquityOption import FinEquityOption
from ...finutils.FinGlobalTypes import FinOptionTypes
from ...market.curves.FinDiscountCurveFlat import FinDiscountCurve
from ...finutils.FinHelperFunctions import labelToString, checkArgumentTypes
from ...finutils.FinDate import FinDate
from ...models.FinModelBlackScholes import bsValue

from scipy.stats import norm
N = norm.cdf

###############################################################################
# TODO: Vectorise pricer
# TODO: NUMBA ??
# TODO: Monte Carlo pricer
###############################################################################


def _f(ss, *args):
    ''' Complex chooser option solve for critical stock price that makes the
    forward starting call and put options have the same price on the chooser
    date. '''

    t = args[0]
    tc = args[1]
    tp = args[2]
    rtc = args[3]
    rtp = args[4]
    kc = args[5]
    kp = args[6]
    v = args[7]
    q = args[8]

    v_call = bsValue(ss, tc-t, kc, rtc, q, v, FinOptionTypes.EUROPEAN_CALL.value)
    v_put = bsValue(ss, tp-t, kp, rtp, q, v, FinOptionTypes.EUROPEAN_PUT.value)

    v = v_call - v_put
    return v

###############################################################################


class FinEquityChooserOption(FinEquityOption):
    ''' A FinEquityChooserOption is an option which allows the holder to
    either enter into a call or a put option on a later expiry date, with both
    strikes potentially different and both expiry dates potentially different.
    This is known as a complex chooser. All the option details are set at trade
    initiation. '''

    def __init__(self,
                 chooseDate: FinDate,
                 callExpiryDate: FinDate,
                 putExpiryDate: FinDate,
                 callStrikePrice: float,
                 putStrikePrice: float):
        ''' Create the FinEquityChooserOption by passing in the chooser date
        and then the put and call expiry dates as well as the corresponding put
        and call strike prices. '''

        checkArgumentTypes(self.__init__, locals())

        if chooseDate > callExpiryDate:
            raise FinError("Expiry date must precede call option expiry date")

        if chooseDate > putExpiryDate:
            raise FinError("Expiry date must precede put option expiry date")

        self._chooseDate = chooseDate
        self._callExpiryDate = callExpiryDate
        self._putExpiryDate = putExpiryDate
        self._callStrike = float(callStrikePrice)
        self._putStrike = float(putStrikePrice)

###############################################################################

    def value(self,
              valueDate: FinDate,
              stockPrice: float,
              discountCurve: FinDiscountCurve,
              dividendCurve: FinDiscountCurve,
              model):
        ''' Value the complex chooser option using an approach by Rubinstein
        (1991). See also Haug page 129 for complex chooser options. '''

        if valueDate > self._chooseDate:
            raise FinError("Value date after choose date.")

        DEBUG_MODE = False

        t = (self._chooseDate - valueDate) / gDaysInYear
        tc = (self._callExpiryDate - valueDate) / gDaysInYear
        tp = (self._putExpiryDate - valueDate) / gDaysInYear

        rt = discountCurve.ccRate(self._chooseDate)
        rtc = discountCurve.ccRate(self._callExpiryDate)
        rtp = discountCurve.ccRate(self._putExpiryDate)

        q = dividendCurve.ccRate(self._chooseDate)

        t = max(t, gSmall)
        tc = max(tc, gSmall)
        tp = max(tp, gSmall)

        v = model._volatility
        v = max(v, gSmall)

        s0 = stockPrice
        xc = self._callStrike
        xp = self._putStrike
        bt = rt - q
        btc = rtc - q
        btp = rtp - q

        argtuple = (t, tc, tp, rtc, rtp, xc, xp, v, q)
        if DEBUG_MODE:
            print("args", argtuple)

        istar = optimize.newton(_f, x0=s0, args=argtuple, tol=1e-8,
                                maxiter=50, fprime2=None)

        if DEBUG_MODE:
            print("istar", istar)

        d1 = (np.log(s0/istar) + (bt + v*v/2)*t) / v / np.sqrt(t)
        d2 = d1 - v * np.sqrt(t)

        if DEBUG_MODE:
            print("d1", d1)
            print("d2", d2)

        y1 = (np.log(s0/xc) + (btc + v*v/2)*tc) / v / np.sqrt(tc)
        y2 = (np.log(s0/xp) + (btp + v*v/2)*tp) / v / np.sqrt(tp)

        if DEBUG_MODE:
            print("y1", y1)
            print("y2", y2)

        rho1 = np.sqrt(t/tc)
        rho2 = np.sqrt(t/tp)

        if DEBUG_MODE:
            print("rho1", rho1)
            print("rho2", rho2)

        w = s0 * np.exp(-q * tc) * M(d1, y1, rho1)
        w = w - xc * np.exp(-rtc * tc) * M(d2, y1 - v * np.sqrt(tc), rho1)
        w = w - s0 * np.exp(-q * tp) * M(-d1, -y2, rho2)
        w = w + xp * np.exp(-rtp * tp) * M(-d2, -y2 + v * np.sqrt(tp), rho2)
        return w

###############################################################################

    def valueMC(self,
                valueDate: FinDate,
                stockPrice: float,
                discountCurve: FinDiscountCurve,
                dividendCurve: FinDiscountCurve,
                model,
                numPaths: int = 10000,
                seed: int = 4242):
        ''' Value the complex chooser option Monte Carlo. '''

        dft = discountCurve.df(self._chooseDate)
        dftc = discountCurve.df(self._callExpiryDate)
        dftp = discountCurve.df(self._putExpiryDate)

        t = (self._chooseDate - valueDate) / gDaysInYear
        tc = (self._callExpiryDate - valueDate) / gDaysInYear
        tp = (self._putExpiryDate - valueDate) / gDaysInYear

        rt = -np.log(dft) / t
        rtc = -np.log(dftc) / tc
        rtp = -np.log(dftp) / tp

        t = max(t, 1e-6)
        tc = max(tc, 1e-6)
        tp = max(tp, 1e-6)

        v = model._volatility
        v = max(v, 1e-6)

        # SHOULD THIS CARE ABOUT TERM STRUCTURE OF Q
        dq = dividendCurve.df(self._chooseDate)
        q = -np.log(dq) / t

#        q = dividendYield
        kc = self._callStrike
        kp = self._putStrike

        np.random.seed(seed)
        sqrtdt = np.sqrt(t)

        # Use Antithetic variables
        g = np.random.normal(0.0, 1.0, size=(1, numPaths))
        s = stockPrice * np.exp((rt - q - v*v / 2.0) * t)
        m = np.exp(g * sqrtdt * v)

        s_1 = s * m
        s_2 = s / m

        v_call_1 = bsValue(s_1, tc-t, kc, rtc, q, v, FinOptionTypes.EUROPEAN_CALL.value)
        v_put_1 = bsValue(s_1, tp-t, kp, rtp, q, v, FinOptionTypes.EUROPEAN_PUT.value)

        v_call_2 = bsValue(s_2, tc-t, kc, rtc, q, v, FinOptionTypes.EUROPEAN_CALL.value)
        v_put_2 = bsValue(s_2, tp-t, kp, rtp, q, v, FinOptionTypes.EUROPEAN_PUT.value)

        payoff_1 = np.maximum(v_call_1, v_put_1)
        payoff_2 = np.maximum(v_call_2, v_put_2)

        payoff = np.mean(payoff_1) + np.mean(payoff_2)
        v = payoff * dft / 2.0
        return v

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("CHOOSER DATE", self._chooseDate)
        s += labelToString("CALL EXPIRY DATE", self._callExpiryDate)
        s += labelToString("CALL STRIKE PRICE", self._callStrike)
        s += labelToString("PUT EXPIRY DATE", self._putExpiryDate)
        s += labelToString("PUT STRIKE PRICE", self._putStrike, "")
        return s

###############################################################################

    def _print(self):
        ''' Simple print function for backward compatibility. '''
        print(self)

###############################################################################
