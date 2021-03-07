###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

import numpy as np
from financepy.models.FinModelSABR import volFunctionSABR
from financepy.models.FinModelSABR import FinModelSABR
from financepy.finutils.FinGlobalTypes import FinOptionTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

alpha = 0.28
beta = 1.0
rho = -0.09
nu = 0.21

f = 0.043
k = 0.050
t = 2.0

###############################################################################


def test_SABR():

    testCases.header("ALPHA", "BETA", "RHO", "VOL")

    for alpha in [0.1, 0.2, 0.3]:
        for beta in [0.5, 1.0, 2.0]:
            for rho in [-0.8, 0.0, 0.8]:
                params = np.array([alpha, beta, rho, nu])
                vol = volFunctionSABR(params, f, k, t)
                testCases.print(alpha, beta, rho, vol)

###############################################################################

def test_SABR_Calibration():

    beta = 0.5
    rho = -0.09
    nu = 0.1
    
    strikeVol = 0.1
    
    f = 0.043
    k = 0.050
    r = 0.03
    texp = 2.0

    callOptionType = FinOptionTypes.EUROPEAN_CALL
    putOptionType = FinOptionTypes.EUROPEAN_PUT

    df = np.exp(-r * texp)

    testCases.header("TEST", "CALIBRATION ERROR")

    # Make SABR equivalent to lognormal (Black) model 
    # (i.e. alpha = 0, beta = 1, rho = 0, nu = 0, shift = 0)
    modelSABR_01 = FinModelSABR(0.0, 1.0, 0.0, 0.0)
    modelSABR_01.setAlphaFromBlackVol(strikeVol, f, k, texp)

    impliedLognormalVol = modelSABR_01.blackVol(f, k, texp)
    impliedATMLognormalVol = modelSABR_01.blackVol(k, k, texp)
    impliedLognormalSmile = impliedLognormalVol - impliedATMLognormalVol

    assert impliedLognormalSmile == 0.0, "In lognormal model, smile should be flat"
    calibrationError = round(strikeVol - impliedLognormalVol, 12)
    testCases.print("LOGNORMAL CASE", calibrationError)
    
    # Volatility: pure SABR dynamics
    modelSABR_02 = FinModelSABR(alpha, beta, rho, nu)
    modelSABR_02.setAlphaFromBlackVol(strikeVol, f, k, texp)

    impliedLognormalVol = modelSABR_02.blackVol(f, k, texp)
    impliedATMLognormalVol = modelSABR_02.blackVol(k, k, texp)
    impliedLognormalSmile = impliedLognormalVol - impliedATMLognormalVol
    calibrationError = round(strikeVol - impliedLognormalVol, 12)
    testCases.print("SABR CASE", calibrationError)

    # Valuation: pure SABR dynamics
    valueCall = modelSABR_02.value(f, k, texp, df, callOptionType)
    valuePut = modelSABR_02.value(f, k, texp, df, putOptionType)
    assert round(valueCall - valuePut, 12) == round(df*(f - k), 12), \
        "The method called 'value()' doesn't comply with Call-Put parity"
        
###############################################################################


test_SABR()
test_SABR_Calibration()

testCases.compareTestCases()
