import numpy as np
from financepy.finutils.FinMath import normcdf_fast, normpdf

from financepy.models.FinModelBlackScholes import bsValue

NC = normcdf_fast
NP = normpdf


class BaronAdesiWhaley:
    """
    American Option Pricing Approximation using Baron Adesi Whaley Model

    S = spot_price
    X = strike_price
    T = expiration_time_in_years
    r = interest_rate_dec_pa
    b = carry_rate_dec_pa
    v = volatility_dec_pa

    """

    def __init__(self, spot_price, strike_price, expiration_time, r, b, vol, phi):
        self.spot_price = spot_price
        self.strike_price = strike_price
        self.expiration_time = expiration_time
        self.rate_of_interest = r
        self.carry_rate = b
        self.volatility = vol
        self.is_call = True if phi==1 else False
        self.E = 0.000001

    def value(self):
        if self.is_call:
            return self._price_american_call()
        else:
            return self._price_american_put()

    def _price_american_call(self):
        # Checked
        if self.carry_rate >= self.rate_of_interest:
            return bsValue(self.spot_price, self.expiration_time, self.strike_price,
                           self.rate_of_interest, self.carry_rate, self.volatility, 1)
        else:
            critical_price = self.get_call_critical_price()
            N = 2 * self.carry_rate / (self.volatility ** 2.0)
            K = 2 * self.rate_of_interest / (
                    (self.volatility ** 2.0) * (1.0 - np.exp(-self.rate_of_interest * self.expiration_time)))
            d1 = (np.log(critical_price / self.strike_price) +
                    (self.carry_rate + (self.volatility ** 2.0) / 2.0) * self.expiration_time) / (
                         self.volatility * np.sqrt(self.expiration_time))
            Q2 = (-1.0 * (N - 1.0) + np.sqrt((N - 1.0) ** 2.0 + 4.0 * K)) / 2.0
            a2 = (critical_price / Q2) * (1.0 - np.exp((self.carry_rate - self.rate_of_interest) * self.expiration_time)
                                          * NC(d1))
            if self.spot_price < critical_price:
                return bsValue(self.spot_price, self.expiration_time, self.strike_price,
                               self.rate_of_interest, self.rate_of_interest - self.carry_rate, self.volatility, 1) + \
                       a2 * ((self.spot_price / critical_price) ** Q2)
            else:
                return self.spot_price - self.strike_price

    def get_call_critical_price(self):
        N = 2.0 * self.carry_rate / (self.volatility ** 2.0)
        M = 2.0 * self.rate_of_interest / (self.volatility ** 2.0)
        q2u = (-(N - 1.0) + np.sqrt(((N - 1.0) * (N - 1.0)) + (4.0 * M))) / 2.0
        su = self.strike_price / (1.0 - (1.0 / q2u))
        h2 = -(self.carry_rate * self.expiration_time + 2.0 * self.volatility * np.sqrt(
            self.expiration_time)) * self.strike_price / (su - self.strike_price)
        Si = self.strike_price + (su - self.strike_price) * (1 - np.exp(h2))
        K = 2 * self.rate_of_interest / (
                (self.volatility ** 2) * (1.0 - np.exp(-self.rate_of_interest * self.expiration_time)))
        d1 = (np.log(Si / self.strike_price) + (self.carry_rate + (self.volatility ** 2.0) / 2.0) * self.expiration_time) / (
                self.volatility * np.sqrt(self.expiration_time))
        Q2 = (-(N - 1.0) + np.sqrt(((N - 1.0)**2.0) + 4.0 * K)) / 2.0
        LHS = Si - self.strike_price
        RHS = bsValue(Si, self.expiration_time, self.strike_price, self.rate_of_interest, self.carry_rate,
                      self.volatility, 1) + \
              (1 - np.exp((self.carry_rate - self.rate_of_interest) * self.expiration_time) * NC(d1)) * Si / Q2
        bi = np.exp((self.carry_rate - self.rate_of_interest) * self.expiration_time) * NC(d1) * (1.0 - 1.0 / Q2)
        bi = bi + (1 - np.exp((self.carry_rate - self.rate_of_interest) * self.expiration_time) * NC(d1) / (
                    self.volatility * np.sqrt(self.expiration_time)))/Q2

        while (np.abs(LHS - RHS) / self.strike_price) < self.E:
            Si = (self.strike_price + RHS - bi * Si) / (1.0 - bi)
            d1 = (np.log(Si / self.strike_price) + (
                        self.carry_rate + (self.volatility ** 2.0) / 2.0) * self.expiration_time) / \
                 (self.volatility * (self.expiration_time ** 0.5))
            LHS = Si - self.strike_price
            RHS = bsValue(Si, self.expiration_time, self.strike_price, self.rate_of_interest, self.carry_rate,
                          self.volatility, 1) + \
                  (1.0 - np.exp((self.carry_rate - self.rate_of_interest) * self.expiration_time) * NC(d1)) * Si / Q2
            bi = np.exp((self.carry_rate - self.rate_of_interest) * self.expiration_time) * NC(d1) * (1.0 - 1.0 / Q2)
            bi = bi + (1 - np.exp((self.carry_rate - self.rate_of_interest) * self.expiration_time) * NP(d1) /
                       (self.volatility * np.sqrt(self.expiration_time)))/Q2
        return Si

    def _price_american_put(self):
        critical_price = self.get_put_critical_price()
        N = 2.0 * self.carry_rate / (self.volatility ** 2.0)
        K = 2.0 * self.carry_rate / (
                    self.volatility ** 2 * (1 - np.exp(-1 * self.rate_of_interest * self.expiration_time)))
        d1 = (np.log(critical_price / self.strike_price) + (self.carry_rate + (self.volatility ** 2) / 2) *
              self.expiration_time) / (self.volatility * (self.expiration_time) ** 0.5)
        Q1 = (-1 * (N - 1) - ((N - 1) ** 2 + 4 * K) ** 0.5) / 2
        a1 = -1 * (critical_price / Q1) * (1 - np.exp((self.carry_rate - self.rate_of_interest)
                                                      * self.expiration_time) * NC(-1 * d1))
        if self.spot_price > critical_price:
            return bsValue(self.spot_price, self.expiration_time, self.strike_price,
                           self.rate_of_interest, 0, self.volatility, -1) + a1 * (
                           self.spot_price / critical_price) * Q1
        else:
            return self.strike_price - self.spot_price

    def get_put_critical_price(self):
        N = 2 * self.carry_rate / self.volatility ** 2
        M = 2 * self.rate_of_interest / self.volatility ** 2
        q1u = (-1 * (N - 1) - ((N - 1) ** 2 + 4 * M) ** 0.5) / 2
        su = self.strike_price / (1 - (1 / q1u))
        h1 = (self.carry_rate * self.expiration_time - 2 * self.volatility * np.sqrt(
            self.expiration_time)) * self.strike_price
        h1 = h1 / (self.strike_price - su)
        Si = su + (self.strike_price - su) * np.exp(h1)
        K = 2 * self.rate_of_interest / (
                self.volatility ** 2 * (1 - np.exp(-self.rate_of_interest * self.expiration_time)))
        d1 = np.log(Si / self.strike_price) + ((self.carry_rate + self.volatility ** 2 / 2) * self.expiration_time) / (
                self.volatility * np.sqrt(self.expiration_time))
        Q1 = (-1 * (N - 1) - ((N - 1) ** 2 + 4 * K) ** 0.5) / 2
        LHS = self.strike_price - Si
        RHS = bsValue(Si, self.expiration_time, self.strike_price, self.rate_of_interest, self.carry_rate,
                      self.volatility, -1) - \
              (1 - np.exp((self.carry_rate - self.rate_of_interest) * self.expiration_time) * NC(-1 * d1)) * Si / Q1
        bi = -1 * np.exp((self.carry_rate - self.rate_of_interest) * self.expiration_time) * NC(-1 * d1) * (1 - 1 / Q1)
        bi = bi - (1 + np.exp((self.carry_rate - self.rate_of_interest) * self.expiration_time)) * NP(-1 * d1) / (
                    self.volatility * np.sqrt(self.expiration_time))
        bi = bi / Q1
        while (np.abs(LHS - RHS) / self.strike_price) < self.E:
            Si = (self.strike_price - RHS + bi * Si) / (1 + bi)
            d1 = (np.log(Si / self.strike_price) + (
                        self.carry_rate + self.volatility ** 2 / 2) * self.expiration_time) / \
                 (self.volatility * (self.expiration_time ** 0.5))
            LHS = self.strike_price - Si
            RHS = bsValue(Si, self.expiration_time, self.strike_price, self.rate_of_interest, self.carry_rate,
                          self.volatility, -1) - \
                  (1 - np.exp((self.carry_rate - self.rate_of_interest) * self.expiration_time) * NC(-1 * d1)) * Si / Q1
            bi = - np.exp((self.carry_rate - self.rate_of_interest) * self.expiration_time) * NC(-1 * d1) * (1 - 1 / Q1)
            bi = bi - (1 + np.exp((self.carry_rate - self.rate_of_interest)) * self.expiration_time) * NP(-1 * d1) / (
                    self.volatility * np.sqrt(self.expiration_time))
            bi = bi / Q1
        return Si


if __name__ == '__main__':
    # spot_price, strike_price, expiration_time, r, b, vol, phi
    baw = BaronAdesiWhaley(100, 100, 0.5, 0.1, 0, 0.35, -1)
    print(baw.value())
