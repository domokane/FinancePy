# -*- coding: utf-8 -*-
"""
Created on Sat Aug  3 14:10:12 2019

@author: Dominic
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Sep 09 14:31:53 2019

@author: Dominic O'Kane
"""
import time
import numpy as np

from financepy.models.FinCIRRateModel import zeroPrice_MC,zeroPrice, FinCIRNumericalScheme
from financepy.finutils.FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__,globalTestCaseMode)

def test_FinCIRRateModel():

    r0 = 0.05
    a = 0.20
    b = 0.05
    sigma = 0.20
    t = 5.0
        
    numPaths = 2000
    dt = 0.05
    seed = 1968
    
    testCases.header("MATURITY","TIME","FORMULA","EULER","LOGNORM","MILSTEIN","KJ","EXACT")

    for t in np.linspace(0,10,21):
    
        start = time.time()      
        p = zeroPrice(r0,a,b,sigma,t)
        p_MC1 = zeroPrice_MC(r0,a,b,sigma,t,dt,numPaths,seed,FinCIRNumericalScheme.EULER.value)
        p_MC2 = zeroPrice_MC(r0,a,b,sigma,t,dt,numPaths,seed,FinCIRNumericalScheme.LOGNORMAL.value)
        p_MC3 = zeroPrice_MC(r0,a,b,sigma,t,dt,numPaths,seed,FinCIRNumericalScheme.MILSTEIN.value)
        p_MC4 = zeroPrice_MC(r0,a,b,sigma,t,dt,numPaths,seed,FinCIRNumericalScheme.KAHLJACKEL.value)
        p_MC5 = zeroPrice_MC(r0,a,b,sigma,t,dt,numPaths,seed,FinCIRNumericalScheme.EXACT.value)    
        end = time.time()
        elapsed = end - start
        testCases.print(t,elapsed,p,p_MC1,p_MC2,p_MC3,p_MC4,p_MC5)
    
    
test_FinCIRRateModel()
testCases.compareTestCases()