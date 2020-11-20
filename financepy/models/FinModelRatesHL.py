##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np
from numba import njit

from ..finutils.FinError import FinError
from ..finutils.FinMath import N
from ..market.curves.FinInterpolator import FinInterpTypes, _uinterpolate
from ..finutils.FinHelperFunctions import labelToString

interp = FinInterpTypes.FLAT_FWD_RATES.value

###############################################################################
# dr = theta(t) dt + sigma * dW
###############################################################################


@njit(fastmath=True, cache=True)
def P_Fast(t, T, Rt, delta, pt, ptd, pT, _sigma):
    ''' Forward discount factor as seen at some time t which may be in the
    future for payment at time T where Rt is the delta-period short rate
    seen at time t and pt is the discount factor to time t, ptd is the one
    period discount factor to time t+dt and pT is the discount factor from
    now until the payment of the 1 dollar of the discount factor. '''

    BtT = (T-t)
    BtDelta = delta
    term1 = np.log(pT/pt) - (BtT/BtDelta) * np.log(ptd/pt)
    term2 = (_sigma**2) * t * BtT * (BtT - BtDelta)/(2.0)

    logAhat = term1 - term2
    BhattT = (BtT/BtDelta) * delta
    p = np.exp(logAhat - BhattT * Rt)
    return p

###############################################################################


class FinModelRatesHL():

    def __init__(self, sigma):
        ''' Construct Ho-Lee model using single parameter of volatility. The
        dynamical equation is dr = theta(t) dt + sigma * dW. Any no-arbitrage
        fitting is done within functions below. '''

        if sigma < 0.0:
            raise FinError("Negative volatility not allowed.")

        self._sigma = sigma

###############################################################################

    def zcb(self, rt1, t1, t2, discountCurve):

        delta = t2 - t1
        dt = 1e-10
        pt1 = discountCurve.df(t1)
        pt1p = discountCurve.df(t1 + dt)
        pt2 = discountCurve.df(t2)
        z = P_Fast(t1, t2, rt1, delta, pt1, pt1p, pt2, self._sigma)
        return z

###############################################################################

    def optionOnZCB(self, texp, tmat,
                    strikePrice, faceAmount,
                    dfTimes, dfValues):
        ''' Price an option on a zero coupon bond using analytical solution of
        Hull-White model. User provides bond face and option strike and expiry
        date and maturity date. '''

        if texp > tmat:
            raise FinError("Option expiry after bond matures.")

        if texp < 0.0:
            raise FinError("Option expiry time negative.")

        ptexp = _uinterpolate(texp, dfTimes, dfValues, interp)
        ptmat = _uinterpolate(tmat, dfTimes, dfValues, interp)

        sigma = self._sigma

        sigmap = sigma * (tmat-texp) * np.sqrt(texp)

        h = np.log((faceAmount*ptmat)/(strikePrice*ptexp))/sigmap+sigmap/2.0

        callValue = faceAmount * ptmat * N(h) - strikePrice * ptexp * N(h-sigmap)
        putValue = strikePrice * ptexp * N(-h+sigmap) - faceAmount * ptmat * N(-h)

        return {'call': callValue, 'put': putValue}

###############################################################################

    def __repr__(self):
        ''' Return string with class details. '''

        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("Sigma", self._sigma)
        return s

###############################################################################
