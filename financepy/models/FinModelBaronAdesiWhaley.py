import numpy as np
from financepy.finutils.FinMath import normcdf_fast, normpdf

from financepy.models.FinModelBlackScholes import bsValue
from scipy import optimize

normcdf = normcdf_fast
NP = normpdf

###############################################################################

def _fcall(si, *args):
    ''' Function to determine ststar for pricing American call options. '''

    t = args[0]
    k = args[1]
    r = args[2]
    q = args[3]
    v = args[4]

    b = r - q
    v2 = v*v

    M = 2.0 * r / v2
    N = 2.0 * b / v2
    K = 1.0 - np.exp(-r * t)

    q2 = (1.0 - N + np.sqrt((N - 1.0)**2 + 4.0 * M/K)) / 2.0
    d1 = (np.log(si / k) + (b + v2 / 2.0) * t) / (v * np.sqrt(t))

    objFn = si - k
    objFn = objFn - bsValue(si, t, k, r, q, v, +1) 
    objFn = objFn - (1.0 - np.exp(-q*t) * normcdf(d1)) * si / q2
    print(objFn)
    return objFn

###############################################################################

def _fput(si, *args):
    ''' Function to determine sstar for pricing American call options. '''

    t = args[0]
    k = args[1]
    r = args[2]
    q = args[3]
    v = args[4]

    b = r - q
    v2 = v*v

    M = 2.0 * r / v2
    N = 2.0 * b / v2
    K = 1.0 - np.exp(-r * t)

    q1 = (1.0 - N - np.sqrt((N - 1.0)**2 + 4.0 * M/K)) / 2.0
    d1 = (np.log(si / k) + (b + v2 / 2.0) * t) / (v * np.sqrt(t))

    objFn = si - k
    objFn = objFn - bsValue(si, t, k, r, q, v, -1) 
    objFn = objFn - (1.0 - np.exp(-q*t) * normcdf(-d1)) * si / q1
    return objFn

###############################################################################


class BaroneAdesiWhaley:
    """
    American Option Pricing Approximation using Baron Adesi Whaley Model

    S = spot_price
    X = strikePrice
    T = timeToExpiry_in_years
    r = interest_rate_dec_pa
    b = carry_rate_dec_pa
    v = volatility_dec_pa

    """

    def __init__(self, 
                 stockPrice: float, 
                 timeToExpiry: float, 
                 strikePrice: float, 
                 riskFreeRate: float, 
                 dividendRate: float, 
                 volatility: float, 
                 phi:int):

        self._stockPrice = stockPrice
        self._strikePrice = strikePrice
        self._timeToExpiry = timeToExpiry
        self._riskFreeRate = riskFreeRate
        self._dividendRate = dividendRate
        self._volatility = volatility
        self._phi = phi

    def value(self):
        if self._phi == 1:
            return self._valueCall()
        else:
            return self._valuePut()

    def _valueCall(self):

        s = self._stockPrice
        t = self._timeToExpiry
        k = self._strikePrice
        r = self._riskFreeRate
        q = self._dividendRate
        b = self._riskFreeRate - self._dividendRate
        v = self._volatility

        if b >= r:
            return bsValue(s, t, k, r, q, v, +1)
        else:

            argtuple = (t, k, r, q, v)

            sstar = optimize.newton(_fcall, x0=s, fprime=None, args=argtuple,
                                    tol=1e-7, maxiter=50, fprime2=None)
                        
            M = 2.0 * r / (v*v)
            N = 2.0 * b / (v*v) 
            K = 1.0 - np.exp(-r * t)
            d1 = (np.log(sstar/k) + (b + v*v/ 2.0) * t) / (v * np.sqrt(t))
            q2 = (-1.0 * (N - 1.0) + np.sqrt((N - 1.0)**2 + 4.0 * M/K)) / 2.0
            A2 = (sstar / q2) * (1.0 - np.exp(-q * t) * normcdf(d1))

            if s < sstar:
                return bsValue(s, t, k, r, q, v, +1) + A2 * ((s/sstar)**q2)
            else:
                return s - k

    def _valuePut(self):
        pass

    # def _price_american_put(self):
    #     critical_price = self.get_put_critical_price()
    #     N = 2.0 * self.carry_rate / (self.volatility ** 2.0)
    #     K = 2.0 * self.carry_rate / (
    #                 self.volatility ** 2 * (1 - np.exp(-1 * self.rate_of_interest * self.timeToExpiry)))
    #     d1 = (np.log(critical_price / self.strikePrice) + (self.carry_rate + (self.volatility ** 2) / 2) *
    #           self.timeToExpiry) / (self.volatility * (self.timeToExpiry) ** 0.5)
    #     Q1 = (-1 * (N - 1) - ((N - 1) ** 2 + 4 * K) ** 0.5) / 2
    #     a1 = -1 * (critical_price / Q1) * (1 - np.exp((self.carry_rate - self.rate_of_interest)
    #                                                   * t) * normcdf(-1 * d1))
    #     if self.spot_price > critical_price:
    #         return bsValue(self.spot_price, self.timeToExpiry, self.strikePrice,
    #                        self.rate_of_interest, 0, self.volatility, -1) + a1 * (
    #                        self.spot_price / critical_price) * Q1
    #     else:
    #         return self.strikePrice - self.spot_price

    # def get_put_critical_price(self):
    #     N = 2 * self.carry_rate / self.volatility ** 2
    #     M = 2 * self.rate_of_interest / self.volatility ** 2
    #     q1u = (-1 * (N - 1) - ((N - 1) ** 2 + 4 * M) ** 0.5) / 2
    #     su = self.strikePrice / (1 - (1 / q1u))
    #     h1 = (self.carry_rate * self.timeToExpiry - 2 * self.volatility * np.sqrt(
    #         self.timeToExpiry)) * self.strikePrice
    #     h1 = h1 / (self.strikePrice - su)
    #     Si = su + (self.strikePrice - su) * np.exp(h1)
    #     K = 2 * self.rate_of_interest / (
    #             self.volatility ** 2 * (1 - np.exp(-self.rate_of_interest * self.timeToExpiry)))
    #     d1 = np.log(Si / self.strikePrice) + ((self.carry_rate + self.volatility ** 2 / 2) * self.timeToExpiry) / (
    #             self.volatility * np.sqrt(self.timeToExpiry))
    #     Q1 = (-1 * (N - 1) - ((N - 1) ** 2 + 4 * K) ** 0.5) / 2
    #     LHS = self.strikePrice - Si
    #     RHS = bsValue(Si, self.timeToExpiry, self.strikePrice, self.rate_of_interest, self.carry_rate,
    #                   self.volatility, -1) - \
    #           (1 - np.exp((self.carry_rate - self.rate_of_interest) * self.timeToExpiry) * NC(-1 * d1)) * Si / Q1
    #     bi = -1 * np.exp((self.carry_rate - self.rate_of_interest) * self.timeToExpiry) * NC(-1 * d1) * (1 - 1 / Q1)
    #     bi = bi - (1 + np.exp((self.carry_rate - self.rate_of_interest) * self.timeToExpiry)) * NP(-1 * d1) / (
    #                 self.volatility * np.sqrt(self.timeToExpiry))
    #     bi = bi / Q1
    #     while (np.abs(LHS - RHS) / self.strikePrice) < self.E:
    #         Si = (self.strikePrice - RHS + bi * Si) / (1 + bi)
    #         d1 = (np.log(Si / self.strikePrice) + (
    #                     self.carry_rate + self.volatility ** 2 / 2) * self.timeToExpiry) / \
    #              (self.volatility * (self.timeToExpiry ** 0.5))
    #         LHS = self.strikePrice - Si
    #         RHS = bsValue(Si, self.timeToExpiry, self.strikePrice, self.rate_of_interest, self.carry_rate,
    #                       self.volatility, -1) - \
    #               (1 - np.exp((self.carry_rate - self.rate_of_interest) * self.timeToExpiry) * NC(-1 * d1)) * Si / Q1
    #         bi = - np.exp((self.carry_rate - self.rate_of_interest) * self.timeToExpiry) * NC(-1 * d1) * (1 - 1 / Q1)
    #         bi = bi - (1 + np.exp((self.carry_rate - self.rate_of_interest)) * self.timeToExpiry) * NP(-1 * d1) / (
    #                 self.volatility * np.sqrt(self.timeToExpiry))
    #         bi = bi / Q1
    #     return Si


if __name__ == '__main__':
    # spot_price, strikePrice, timeToExpiry, r, b, vol, phi

    # Checking against table 3-1 in Haug 
    k = 100.0
    r = 0.10
    q = 0.10
    
    for t in [0.1, 0.5]:
        for v in [0.15, 0.25, 0.35]:
            for s in [90.0, 100.0, 110.0]:
                
                baw = BaroneAdesiWhaley(s, t, k, r, q, v, 1)
                print("%9.5f %9.5f %9.5f %9.5f"% (s, t, v, baw.value()))

