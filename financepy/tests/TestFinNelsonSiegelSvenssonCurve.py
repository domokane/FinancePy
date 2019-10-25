# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 09:26:27 2016

@author: Dominic O'Kane
"""
import numpy as np
import matplotlib.pyplot as plt

from financepy.finutils.FinMath import scale
from financepy.market.curves.FinNelsonSiegelSvenssonCurve import FinNelsonSiegelSvenssonCurve

from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__,globalTestCaseMode)

################################################################################

def test_FinNelsonSiegelSvenssonCurve():

    tau1 = 2.0
    tau2 = 0.5
    times = np.linspace(0.0,10.0,5)

    curve1 = FinNelsonSiegelSvenssonCurve(1.,0.,0.,0.,tau1,tau2)
    factor1loading = curve1.zero(times)
    curve2 = FinNelsonSiegelSvenssonCurve(0.,1.,0.,0.,tau1,tau2)
    factor2loading = curve2.zero(times)
    curve3 = FinNelsonSiegelSvenssonCurve(0.,0.,1.,0.,tau1,tau2)
    factor3loading = curve3.zero(times)
    curve4 = FinNelsonSiegelSvenssonCurve(0.,0.,0.,1.,tau1,tau2)
    factor4loading = curve4.zero(times)

    testCases.header("FACTOR LOADING ON ZERO RATES")
    testCases.print(factor1loading)
    testCases.print(factor2loading)
    testCases.print(factor3loading)
    testCases.print(factor4loading)

#    plt.figure(figsize = (6,4))
#    plt.plot(times,scaleVector(factor1loading,1),label='beta1');
#    plt.plot(times,scaleVector(factor2loading,1),label='beta2');
#    plt.plot(times,scaleVector(factor3loading,1),label='beta3');
#    plt.ylim((0,1.05))
#
#    plt.title('Factor Loadings in Nelson-Siegel Model');
#    plt.xlabel('Time (years)');
#    plt.ylabel('Loading');
#    plt.legend(loc='best')

################################################################################

    testCases.header("BETA1","BETA2","BETA3","BETA4","ZEROS")

    beta1 = 0.03
    beta2 = -0.02
    beta3 = -0.02
    beta4 = 0.08
    curve1 = FinNelsonSiegelSvenssonCurve(beta1,beta2,beta3,beta4,tau1,tau2)
    zeroRates1 = curve1.zero(times);
    testCases.print(beta1,beta2,beta3,beta4,zeroRates1)

    beta1 = 0.04
    beta2 = -0.02
    beta3 = -0.02
    beta4 = 0.08
    curve2 = FinNelsonSiegelSvenssonCurve(beta1,beta2,beta3,beta4,tau1,tau2)
    zeroRates2 = curve2.zero(times);
    testCases.print(beta1,beta2,beta3,beta4,zeroRates2)

    beta1 = 0.05
    beta2 = -0.02
    beta3 = -0.02
    beta4 = 0.08
    curve3 = FinNelsonSiegelSvenssonCurve(beta1,beta2,beta3,beta4,tau1,tau2)
    zeroRates3 = curve3.zero(times);
    testCases.print(beta1,beta2,beta3,beta4,zeroRates3)

    beta1 = 0.06
    beta2 = -0.02
    beta3 = -0.02
    beta4 = 0.08
    curve4 = FinNelsonSiegelSvenssonCurve(beta1,beta2,beta3,beta4,tau1,tau2)
    zeroRates4 = curve4.zero(times);
    testCases.print(beta1,beta2,beta3,beta4,zeroRates4)

    beta1 = 0.07
    beta2 = -0.02
    beta3 = -0.02
    beta4 = 0.08
    curve5 = FinNelsonSiegelSvenssonCurve(beta1,beta2,beta3,beta4,tau1,tau2)
    zeroRates5 = curve5.zero(times);
    testCases.print(beta1,beta2,beta3,beta4,zeroRates5)

    if 1==0:
        plt.figure(figsize = (6,4))
        plt.plot(times,scale(zeroRates1,100),label='beta1=3%');
        plt.plot(times,scale(zeroRates2,100),label='beta1=4%');
        plt.plot(times,scale(zeroRates3,100),label='beta1=5%');
        plt.plot(times,scale(zeroRates4,100),label='beta1=6%');
        plt.plot(times,scale(zeroRates5,100),label='beta1=7%');
        plt.ylim((0,7.5))
    
        plt.title('Nelson-Siegel Zero Rate Curves');
        plt.xlabel('Time (years)');
        plt.ylabel('Zero Rate (%)');
        plt.legend(loc='lower right',frameon=False)

################################################################################

test_FinNelsonSiegelSvenssonCurve()
testCases.compareTestCases()