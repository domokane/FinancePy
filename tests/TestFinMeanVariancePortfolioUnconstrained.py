"""
This script serves to show the link between mean-variance portfolios

1) Compares the mean-variance weights from Markowitz MV to Hansen-Richard MV
2) Computes the mean-variance weights based on a Factor model and factor-
mimicking porfolios to mean-variance weights based on Markowitz MV with a
reduced covariance matrix.

"""
import sys

import numpy as np
import pandas as pd

import financepy.portfolio as fpp
from FinTestCases import FinTestCases, globalTestCaseMode

sys.path.append("..//..")

testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################
# get the data

# The data is retrieved from the excellent Data page by Prof. Ken French
# https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html

###############################################################################


def test_MeanVarianceOptimisation():

    ff_factors = pd.read_pickle('data/ff_factors.pkl')
    ff_assets = pd.read_pickle('data/ff_assets.pkl')

    num_times = len(ff_factors)

    # descriptive statistics
    print('(Monthly) Expected returns Fama-French Factors\n',
          ff_factors.mean())
    print('(Monthly) Expected returns Fama-French assets\n',
          ff_assets.mean())

    print('(Monthly) volatility Fama-French Factors\n', ff_factors.std())
    print('(Monthly) volatility Fama-French assets\n', ff_assets.std())

    ###########################################################################
    #
    # Compute expected returns and covariance matrix
    #
    ###########################################################################

    # compute R
    # Note that Markowitz worked with (gross) arithmetic returns
    R_norf = ff_assets + 1
    # drop the risk-free rate
    R_norf.drop('rf', axis=1, inplace=True)

    # compute estimates of average returns per asset class
    # 1/T \sum R^{i}_t
    mu_norf = R_norf.mean(0)
    # compute covariance matrix of R
    # this is an unbiased version [e.g 1 / (T-1)]
    cov_norf = R_norf.cov()

    ###########################################################################
    #
    # Compute Classical Mean-Variance portfolio Markowitz
    #
    ###########################################################################

    # mean-variance instance based on Markowitz quadratic optimization
    mvm = fpp.MeanVarianceMarkowitz(mu_norf,
                                    cov_norf * (num_times - 1) / num_times)

    mus = np.arange(-100, 100.1) / 1000


    # computes the expected return for the minimum variance attainable
    mvm_min_var = mvm.rate_min_var()

    # return the variance for each desired expected return
    mvm_vars = mvm.mv_frontier(mus + mvm_min_var)

    ###########################################################################
    #
    # Compute Modern Mean-Variance portfolio Hansen-Richard
    #
    ###########################################################################

    # mean-variance instance on Hansen-Richard orthogonal decomposition
    mvhr = fpp.MeanVarianceHansenRichard(mu_norf,
                                         cov_norf * (num_times - 1)
                                         / num_times)

    # computes the expected return for the minimum variance attainable
    mvhr_min_var = mvhr.rate_min_var()

    # return the variance for each desired expected return
    mvhr_vars = mvhr.mv_frontier(mus + mvhr_min_var)

    testCases.header("MVM_VARS", "MV_HR_VARS")
    testCases.print(mvm_vars, mvhr_vars)

###############################################################################

test_MeanVarianceOptimisation()
testCases.compareTestCases()
