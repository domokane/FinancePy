##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import numpy as np

from ..utils.math import N, M, phi3
from ..utils.math import norminvcdf as NormSInv
from ..utils.error import FinError

###############################################################################


class LHPlusModel():
    """ Large Homogenous Portfolio model with extra asset. Used for
    approximating full Gaussian copula. """

    def __init__(self, p, r, h, beta, p0, r0, h0, beta_0):
        self._p = p
        self._r = r
        self._h = h
        self._beta = beta
        self._p0 = p0
        self._r0 = r0
        self._h0 = h0
        self._beta_0 = beta_0

###############################################################################

    def prob_loss_gt_k(self, k):
        """ Returns P(L>K) where L is the portfolio loss given by model. """

        if k < (1.0 - self._r0) * self._h0:
            raise FinError("Function does not work when K<(1-R0)H0")

        c = NormSInv(self._p)
        c0 = NormSInv(self._p0)
        arga = K / (1.0 - self._r) / self._h
        inva = 0.0

        if arga < 0.00001:
            inva = -9999999999999
        elif arga >= 1:
            inva = 9999999999999.0
        else:
            inva = NormSInv(arga)

        RtOneMinusBeta2 = np.sqrt(1.0 - self._beta * self._beta)
        a = (1.0 / self._beta) * (c - RtOneMinusBeta2 * inva)
        argb = (K - (1.0 - self._r0) * self._h0) / (1.0 - self._r) / self._h
        invb = 0.0

        if argb <= 0:
            invb = -9999999999999.0
        elif argb > 1.0:
            invb = 99999999999999.0
        else:
            invb = NormSInv(argb)

        b = (1.0 / self._beta) * \
            ((c - RtOneMinusBeta2 * invb))

        prob_loss_gt_k = N(a) + M(c0, b, self._beta_0) - M(c0, a, self._beta_0)

        return prob_loss_gt_k

###############################################################################

    def exp_min_lk_integral(self, k, dk):

        k0 = 0.0
        num_steps = int(k / dk)
        dk = k / num_steps
        cdf0 = 1.0
        cdf1 = 0.0
        exp_min_lk = 0.0
        check_sum = 0.0

        for _ in range(0, num_steps):
            k0 += dk
            cdf1 = self.prob_loss_gt_k(k0)
            pdf = cdf0 - cdf1
            cdf0 = cdf1
            check_sum += pdf
            exp_min_lk += pdf * k0

        check_sum += cdf1
        exp_min_lk += cdf1 * K

        return exp_min_lk

###############################################################################

    def exp_min_lk(self, k):

        if abs(k) < 1e-6:
            return k

        c = NormSInv(self._p)
        c0 = NormSInv(self._p0)
        arga = k / (1.0 - self._r) / self._h
        argb = (k - (1.0 - self._r0) * self._h0) / (1.0 - self._r) / self._h

        if argb < 0:
            raise FinError("Tranche too thin for LHPlus")

        inva = NormSInv(arga)
        invb = NormSInv(argb)
        rt_one_minus_beta2 = np.sqrt(1.0 - self._beta * self._beta)
        a = (c - rt_one_minus_beta2 * inva) / self._beta
        b = (c - rt_one_minus_beta2 * invb) / self._beta

        r12 = self._beta
        r13 = self._beta_0
        r23 = self._beta * self._beta_0
        el1 = self._p * self._h * (1.0 - self._r)
        el2 = self._p0 * self._h0 * (1.0 - self._r0)
        el3 = -K * (M(c0, a, self._beta_0) - N(a))
        el4 = - ((1.0 - self._r0) * self._h0 - k) * M(c0, b, self._beta_0)
        term1 = M(c, a, self._beta) + phi3(b, c, c0, r12, r13, r23) \
            - phi3(a, c, c0, r12, r13, r23)
        el5 = - (1.0 - self._r) * self._h * term1
        elk1k2 = el1 + el2 + el3 + el4 + el5
        return elk1k2

###############################################################################

    def exp_min_lk2(self, k):

        if abs(k) < 1e-6:
            return k

        c = NormSInv(self._p)
        c0 = NormSInv(self._p0)
        arga = k / (1.0 - self._r) / self._h
        argb = (k - (1.0 - self._r0) * self._h0) / (1.0 - self._r) / self._h

        if argb < 0.0:
            raise FinError("Tranche too thin for LHPlus")

        inva = NormSInv(arga)
        invb = NormSInv(argb)

        RtOneMinusBeta2 = np.sqrt(1.0 - self._beta * self._beta)

        a = (1.0 / self._beta) * (c - RtOneMinusBeta2 * inva)
        b = (1.0 / self._beta) * (c - RtOneMinusBeta2 * invb)

        r12 = self._beta
        r13 = self._beta_0
        r23 = self._beta * self._beta_0

        el1 = self._p * self._h * (1.0 - self._r) + \
            self._p0 * self._h0 * (1.0 - self._r0)
        el1 = el1 - k * (M(c0, a, self._beta_0) - N(a))
        el1 = el1 - ((1.0 - self._r0) * self._h0 - k) * M(c0, b, self._beta_0)

        term = M(c, a, self._beta) + phi3(b, c, c0, r12, r13, r23) \
            - phi3(a, c, c0, r12, r13, r23)

        el1 = el1 - (1.0 - self._r) * self._h * term
        return el1

###############################################################################

    def tranche_survival_prob(self, k1, k2):

        if k2 == k1:
            raise FinError("tranche_survival_prob: Same strikes")

        dk = 0.00001
        el_k2 = self.exp_min_lk_integral(k2, dk)
        el_k1 = self.exp_min_lk_integral(k1, dk)
        q = 1.0 - (el_k2 - el_k1) / (k2 - k1)
        return q

###############################################################################
