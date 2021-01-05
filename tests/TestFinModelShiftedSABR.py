###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
sys.path.append("..")

import numpy as np
from financepy.models.FinModelSABRShifted import FinModelSABRShifted
from financepy.finutils.FinGlobalTypes import FinOptionTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################

beta = 0.5
rho = -0.09
nu = 0.1
shift = 0.02

strikeVol = 0.1

forward = 0.043
strike = 0.050
riskFreeIR = 0.03
timeToExpiry = 1.0

df = np.exp(-riskFreeIR * timeToExpiry)

callOptionType = FinOptionTypes.EUROPEAN_CALL
putOptionType = FinOptionTypes.EUROPEAN_PUT

###############################################################################


def test_ShiftedSABR():

    testCases.header("TEST", "CALIBRATION ERROR")

    # Volatility: SABR equivalent to lognormal (Black) model (i.e. beta = 1, rho = 0, nu = 0, shift = 0)
    modelSABR_01 = FinModelSABRShifted(1.0, 0, 0.0, 0.0)
    modelSABR_01.calibrateAlpha(strikeVol, 'LOGNORMAL', forward, strike, timeToExpiry)

    impliedLognormalVol = modelSABR_01.implyLognormalVol(forward, strike, timeToExpiry)
    impliedATMLognormalVol = modelSABR_01.implyLognormalVol(strike, strike, timeToExpiry)
    impliedLognormalSmile = impliedLognormalVol - impliedATMLognormalVol
    assert impliedLognormalSmile == 0.0, "In lognormal model, smile should be flat"
    calibrationError = round(strikeVol - impliedLognormalVol, 12)
    testCases.print("LOGNORMAL CASE", calibrationError)

    # Volatility: pure SABR dynamics
    modelSABR_02 = FinModelSABRShifted(beta, rho, nu, shift)
    modelSABR_02.calibrateAlpha(strikeVol, 'LOGNORMAL', forward, strike, timeToExpiry)

    impliedLognormalVol = modelSABR_02.implyLognormalVol(forward, strike, timeToExpiry)
    impliedATMLognormalVol = modelSABR_02.implyLognormalVol(strike, strike, timeToExpiry)
    impliedLognormalSmile = impliedLognormalVol - impliedATMLognormalVol
    calibrationError = round(strikeVol - impliedLognormalVol, 12)
    testCases.print("SABR CASE", calibrationError)

    # Valuation: pure SABR dynamics
    valueCall = modelSABR_02.value(forward, strike, riskFreeIR, timeToExpiry, callOptionType)
    valuePut = modelSABR_02.value(forward, strike, riskFreeIR, timeToExpiry, putOptionType)
    assert round((1 / df) * (valueCall - valuePut), 12) == round(forward - strike, 12), \
        "The method called 'value()' doesn't comply with Call-Put parity"

    # TODO: adding Call-Put parity test for all sensitivities

###############################################################################


test_ShiftedSABR()
testCases.compareTestCases()