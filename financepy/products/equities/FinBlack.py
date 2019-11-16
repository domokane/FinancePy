# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 15:28:41 2019

@author: Dominic O'Kane
"""
# TODO Fix this

import numpy as np
from scipy.stats import norm

from ...finutils.FinMath import N


class BlackModel():

    ##########################################################################

    def value(self,
              forwardRate,
              strikeRate,
              timeToExpiry,
              sigma,
              callOrPut):

        if strikeRate < 0:
            raise Exception("Error: Negative strike")

        if strikeRate == 0.0:
            strike = 1e-16

        sqrtT = np.sqrt(timeToExpiry)

        d1 = np.log(forwardRate / strike) + sigma * sigma * timeToExpiry / 2
        d1 = d1 / (sigma * sqrtT)
        d2 = d1 - sigma * sqrtT

        returnValue = 0.0

        if callOrPut == "C":
            returnValue = forwardRate * norm.cdf(d1) - strike * N(d2)
        elif callOrPut == "P":
            returnValue = strike * norm.cdf - forwardRate * N(-d1)
        else:
            raise Exception("Option type must be a Call(C) or a Put(P)")

        return returnValue

###############################################################################

#    public static double BlackImpliedVol(double forward,
#                                         double timeToExpiry,
#                                         double strike,
#                                         ref string callOrPut,
#                                         double optionPrice)
#    {
#
#      BlackVolObjectiveFunction objFnBlackModel = new BlackVolObjectiveFunction();
#      objFnBlackModel.m_forward = forward;
#      objFnBlackModel.m_strike = strike;
#      objFnBlackModel.m_timeToExpiry = timeToExpiry;
#      objFnBlackModel.m_callOrPut = callOrPut;
#      objFnBlackModel.m_optionPrice = optionPrice;
#
#      // optimisation parameters
#      double xmin = 0.00000001;
#      double xmax = 0.9999;
#      double xacc = 0.0000000001;
#      int maxIter = 40;
#
#      ObjectiveFunction objFn = (ObjectiveFunction) objFnBlackModel;
#      Solvers.generalOneDimensionalSolver(ref objFn, xmin, xmax, xacc, xacc, maxIter);
#
#      return objFnBlackModel.m_volatility;
#    }
#
