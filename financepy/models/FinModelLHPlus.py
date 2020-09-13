##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ..finutils.FinMath import N, M, phi3
from ..finutils.FinMath import norminvcdf as NormSInv
from ..finutils.FinError import FinError

###############################################################################


class LHPlusModel():
    ''' Large Homogenous Portfolio model with extra asset. Used for
    approximating full Gaussian copula. '''

    def __init__(self, P, R, H, beta, P0, R0, H0, beta0):
        self._P = P
        self._R = R
        self._H = H
        self._beta = beta
        self._P0 = P0
        self._R0 = R0
        self._H0 = H0
        self._beta0 = beta0

###############################################################################

    def probLossGreaterThanK(self, K):
        ''' Returns P(L>K) where L is the portfolio loss given by model. '''
        if K < (1.0 - self._R0) * self._H0:
            raise FinError("Function does not work when K<(1-R0)H0")

        c = NormSInv(self._P)
        c0 = NormSInv(self._P0)
        arga = K / (1.0 - self._R) / self._H
        inva = 0.0

        if arga < 0.00001:
            inva = -9999999999999
        elif arga >= 1:
            inva = 9999999999999.0
        else:
            inva = NormSInv(arga)

        RtOneMinusBeta2 = np.sqrt(1.0 - self._beta * self._beta)
        a = (1.0 / self._beta) * (c - RtOneMinusBeta2 * inva)
        argb = (K - (1.0 - self._R0) * self._H0) / (1.0 - self._R) / self._H
        invb = 0.0

        if argb <= 0:
            invb = -9999999999999.0
        elif argb > 1.0:
            invb = 99999999999999.0
        else:
            invb = NormSInv(argb)

        b = (1.0 / self._beta) * \
            ((c - RtOneMinusBeta2 * invb))

        probLgtK = N(a) + M(c0, b, self._beta0) - M(c0, a, self._beta0)

        return probLgtK

###############################################################################

    def expMinLKIntegral(self, K, dK):

        k0 = 0.0
        numSteps = int(K / dK)
        dK = K / numSteps
        cdf0 = 1.0
        cdf1 = 0.0
        expMinLK = 0.0
        checkSum = 0.0

        for _ in range(0, numSteps):
            k0 += dK
            cdf1 = self.probLossGreaterThanK(k0)
            pdf = cdf0 - cdf1
            cdf0 = cdf1
            checkSum += pdf
            expMinLK += pdf * k0

        checkSum += cdf1
        expMinLK += cdf1 * K

        return expMinLK

###############################################################################

    def expMinLK(self, K):

        if abs(K) < 1e-6:
            return K

        c = NormSInv(self._P)
        c0 = NormSInv(self._P0)
        arga = K / (1.0 - self._R) / self._H
        argb = (K - (1.0 - self._R0) * self._H0) / (1.0 - self._R) / self._H

        if argb < 0:
            raise FinError("Tranche too thin for LHPlus")

        inva = NormSInv(arga)
        invb = NormSInv(argb)

        RtOneMinusBeta2 = np.sqrt(1.0 - self._beta * self._beta)

        a = 1.0 / self._beta * (c - RtOneMinusBeta2 * inva)
        b = 1.0 / self._beta * (c - RtOneMinusBeta2 * invb)

        r12 = self._beta
        r13 = self._beta0
        r23 = self._beta * self._beta0

        el1 = self._P * self._H * (1.0 - self._R)
        el2 = self._P0 * self._H0 * (1.0 - self._R0)
        el3 = -K * (M(c0, a, self._beta0) - N(a))
        el4 = - ((1.0 - self._R0) * self._H0 - K) * M(c0, b, self._beta0)
        term1 = M(c, a, self._beta) + phi3(b, c, c0, r12, r13, r23) \
            - phi3(a, c, c0, r12, r13, r23)
        el5 = - (1.0 - self._R) * self._H * term1

        elk1k2 = el1 + el2 + el3 + el4 + el5
        return elk1k2

###############################################################################

    def expMinLK2(self, K):

        if abs(K) < 1e-6:
            return K

        c = NormSInv(self._P)
        c0 = NormSInv(self._P0)

        arga = K / (1.0 - self._R) / self._H
        argb = (K - (1.0 - self._R0) * self._H0) / (1.0 - self._R) / self._H

        if argb < 0.0:
            raise FinError("Tranche too thin for LHPlus")

        inva = NormSInv(arga)
        invb = NormSInv(argb)

        RtOneMinusBeta2 = np.sqrt(1.0 - self._beta * self._beta)

        a = (1.0 / self._beta) * (c - RtOneMinusBeta2 * inva)
        b = (1.0 / self._beta) * (c - RtOneMinusBeta2 * invb)

        r12 = self._beta
        r13 = self._beta0
        r23 = self._beta * self._beta0

        el1 = self._P * self._H * (1.0 - self._R) + \
            self._P0 * self._H0 * (1.0 - self._R0)
        el1 = el1 - K * (M(c0, a, self._beta0) - N(a))
        el1 = el1 - ((1.0 - self._R0) * self._H0 - K) * M(c0, b, self._beta0)

        term = M(c, a, self._beta) + phi3(b, c, c0, r12, r13, r23) \
            - phi3(a, c, c0, r12, r13, r23)

        el1 = el1 - (1.0 - self._R) * self._H * term
        return el1

###############################################################################

    def trancheSurvivalProbability(self, k1, k2):

        if k2 == k1:
            raise FinError("trancheSurvivalProbability: Same strikes")

        dk = 0.00001
        elK2 = self.expMinLKIntegral(k2, dk)
        elK1 = self.expMinLKIntegral(k1, dk)
        q = 1.0 - (elK2 - elK1) / (k2 - k1)
        return q

###############################################################################
