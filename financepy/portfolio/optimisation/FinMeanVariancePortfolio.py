###############################################################################
# Copyright (C) 2018, 2019, 2020 Jeroen Kerkhof
###############################################################################


###############################################################################

"""
Mean-Variance classes based on Mean-Variance Analysis in Cochrane(2005)
"""

import numpy as np
import pandas as pd
import numpy.linalg as la

###############################################################################


class MeanVarianceMarkowitz:
    """
    Mean-Variance class based on Markowitz quadratic programming
    """
    def __init__(self, mu, sigma, rf=None):
        """
        Args:
            mu: expected gross returns
            cov: covariance matrix of gross returns
        """
        self.num_assets = len(mu)
        self.mu = mu
        self.sigma = sigma
        self.rf = rf
        # vector of ones
        iota_n = np.ones((len(mu),))
        self.sigma_inv = la.inv(sigma)
        # helper variables [see p. 82 Cochrane]
        self.A = self.mu @ self.sigma_inv @ self.mu
        self.B = self.mu @ self.sigma_inv @ iota_n
        self.C = iota_n @ self.sigma_inv @ iota_n

        # weight are a combination of weights of two portfolios
        # e.g. Mutual Fund Separation Theorem [e.g. Campbell(2017, 2.57)]
        # Though you would need to rewrite it a bit
        # as w_i and w_mu are not proper weights as they do not add to 1.
        self.base_w_i = self.sigma_inv @ iota_n
        self.base_w_mu = self.sigma_inv @ mu

    def min_var(self):
        """
        the lowest possible variance
        """
        return self.B / self.C

    def mv_frontier(self, mus):
        """
        computes the mean-variance frontier
        using the analytic solution to the quadratric programming
        problem.
        """
        # Cochrane (2005) (5.7)
        mv_vars = self.A * np.ones(len(mus))
        for i, m in enumerate(mus):
            mv_vars[i] += (self.C*m**2 - 2*self.B*m)

        mv_vars /= (self.A*self.C - self.B**2)

        return mv_vars

    def mv_weights(self, mu):
        """
        compute the mean-variance weights for a given expected return
        """
        w = (self.base_w_mu * (self.C * mu - self.B) +
             self.base_w_i * (self.A - self.B*mu)) / \
             (self.A * self.C -self.B**2)

        return w

    def tangency_portfolio(self, rf):
        """
        computes the tangency porfolio for a given risk-free rate
        """
        # Campbell(2017, 2.62)
        excess_r = (self.mu - rf*np.ones(len(self.mu)))
        w = self.sigma_inv @ excess_r
        # rescale weight to add to 1
        w /= w.sum()

        return w

###############################################################################


class MeanVarianceHansenRichard:
    """
    This class provides the Mean-Variance calculations for the Hansen-Richard
    version of the Mean-variance frontier.

    It assumes there is no risk-free asset
    """

    def __init__(self, mu, cov):
        """
        Args:
            mu: expected gross returns
            cov: covariance matrix of gross returns
        """
        self.num_assets = len(mu)
        # TODO add shape check on cov
        # mean of the returns
        self.mu = mu
        # covariance matrix of the returns
        self.cov = cov
        # second moments
        ERR = cov + np.outer(mu, mu)
        # inverse of second moments
        self.ERR_inv = la.inv(ERR)

        # Compute R^{*}
        # p.91 Cochrane(2005) (13)
        iota_n = np.ones(self.num_assets)
        r_star_weights = iota_n @ self.ERR_inv
        r_star_weights /= iota_n @ self.ERR_inv @ iota_n
        self.r_star_weights = pd.Series(r_star_weights,
                                        index=self.mu.index)

        # compute the expected return of R^{*}
        self.e_r_star = r_star_weights @ self.mu

        # Compute R^{e*}
        # We use R^{e} = R - R^{*}
        self.mu_e = self.mu - self.e_r_star

        # Compute second moments of excess returns
        EReRe = ERR - r_star_weights @ ERR @ r_star_weights
        # project the unit payoff on the excess return space
        # compute the inverse of the second moment. Note rank = n-1
        self.EReRe_inv = la.pinv(EReRe)
        self.re_star_weights = pd.Series(self.mu_e @ self.EReRe_inv,
                                         index=self.mu.index)

        # compute E[R^{e*}]
        self.e_re_star = self.re_star_weights @ self.mu_e

        # compute E[R^{*2}]
        self.e_r_star2 = r_star_weights @ ERR @ r_star_weights

        self.weights_adjustment = (np.eye(self.num_assets) -
                                   np.outer(self.r_star_weights, iota_n)) @ \
                                   self.re_star_weights

    def r_star_weights(self):
        """
        weights of R^{*}, the return with the smallest second moment
        """
        return self.r_star_weights

    def re_star_weights(self):
        """
        weights of R^{e*}, the projection of unit payoff on excess returns
        """
        return self.re_star_weights

    def e_r_star(self):
        """
        the expected return of R^{*}
        """
        return self.e_r_star

    def _weight_zero_beta(self):
        """
        the weight for the zero-beta return

        TODO check unbiased vs biased variance estimate
        """
        return (self.e_r_star2 - self.e_r_star**2) / \
            (self.e_r_star * self.e_re_star)

    def _weight_min_var(self):
        """
        the weight for the minimum variance portfolio
        """
        return self.e_r_star / (1.0 - self.e_re_star)

    def _weight_constant_mimicking_portfolio(self):
        """
        the weight for the constant mimicking portfolio
        """
        return self.e_r_star2 / self.e_r_star

    def rate_min_var(self):
        """
        the minimum variance
        """
        return self.e_r_star + self._weight_min_var() * self.e_re_star

    def rate_zero_beta(self):
        """
        the zero beta rate Cochrane (2005) p.112
        Note: the zero beta rate equals the weight of the constant
        mimicking portfolio.
        """
        return self._weight_constant_mimicking_portfolio()

    def rate_constant_mimicking_portfolio(self):
        """
        the rate on the constant mimicking portfolio
        """
        return self.e_r_star + self._weight_constant_mimicking_portfolio() \
            * self.e_re_star

    def mv_frontier(self, mus):
        """
        computes the mean-variance frontier

        it computes the mean-variance variances for given expected returns

        Args:
            mus (array-like): expected returns on the MV frontier

        Returns:
            mv_vars (array-like): associated variances of the MV frontier
        """
        mv_vars = np.zeros(len(mus))
        # p.114 Cochrane (2005)
        for i, mu in enumerate(mus):
            w = mu / self.e_re_star + self._weight_min_var()
            mv_vars[i] = self.e_r_star2 + w**2 * self.e_re_star \
                - self.e_r_star**2 - 2*w*self.e_r_star*self.e_re_star \
                - w**2*self.e_re_star**2

        return mv_vars

    def mv_weights(self, mu):
        """
        compute the mean-variance weights for a given expected return
        """
        # compute the weight on the HR decomposition
        w_mu = (mu - self.e_r_star) / self.e_re_star
        # calculate the weights
        weights = self.r_star_weights + w_mu * self.weights_adjustment

        return weights

    def hj_bounds(self, range_em=np.arange(90, 110) / 100):
        """
        compute the Hansen-Jaganathan bounds
        see e.g. Ferson(2019) (10.17)
        """
        # create a vector of ones
        iota_n = np.ones((len(self.mu),))
        # loop over all given values for E[m]
        v = np.zeros(len(range_em))
        for i, em in enumerate(range_em):
            one_minus_emer = (iota_n - em*self.mu)
            v[i] = one_minus_emer.T @ la.inv(self.cov) @ one_minus_emer

        return v
