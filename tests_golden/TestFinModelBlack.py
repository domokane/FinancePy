###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
# Guillaume Lefieux
###############################################################################

import sys
import numpy as np
sys.path.append("..")

from financepy.models.FinModelBlack import FinModelBlack
from financepy.finutils.FinGlobalTypes import FinOptionTypes

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)


def test_Black():

    forward = 0.034
    strike = 0.050
    riskFreeIR = 0.00
    timeToExpiry = 2.0
    volatility = 0.20

    testCases.header("ITEM", "CALL", "PUT")

    callOptionType = FinOptionTypes.EUROPEAN_CALL
    putOptionType = FinOptionTypes.EUROPEAN_PUT

    df = np.exp(-riskFreeIR * timeToExpiry)
    model = FinModelBlack(volatility)

    dp = 12 # Precision
    
    try:

        #######################################################################

        valueCall = model.value(forward, strike, timeToExpiry, df, callOptionType)
        valuePut = model.value(forward, strike, timeToExpiry, df, putOptionType)

        assert round((valueCall - valuePut), dp) == round(df*(forward - strike), dp), \
            "The method called 'value()' doesn't comply with Call-Put parity"

        testCases.print("VALUE", valueCall, valuePut)

        #######################################################################

        deltaCall = model.delta(forward, strike, timeToExpiry, df, callOptionType)
        deltaPut = model.delta(forward, strike, timeToExpiry, df, putOptionType)

        assert round((1/df) * (deltaCall - deltaPut), dp) == 1.0, \
            "The method called 'delta()' doesn't comply with Call-put parity"

        testCases.print("DELTA", deltaCall, deltaPut)

        #######################################################################

        gammaCall = model.gamma(forward, strike,timeToExpiry, df, callOptionType)
        gammaPut = model.gamma(forward, strike, timeToExpiry, df, putOptionType)

        assert round(gammaCall - gammaPut, dp) == 0.0, \
            "The method called 'gamma()' doesn't comply with Call-Put parity"

        testCases.print("GAMMA", gammaCall, gammaPut)

        #######################################################################

        thetaCall = model.theta(forward, strike, timeToExpiry, df, callOptionType)
        thetaPut = model.theta(forward, strike, timeToExpiry, df, putOptionType)

        assert round((thetaCall - thetaPut), dp) == round((riskFreeIR * timeToExpiry) * (forward - strike) * df, dp), \
            "The method called 'theta()' doesn't comply with Call-Put parity"

        testCases.print("THETA", thetaCall, thetaPut)

        #######################################################################

        vegaCall = model.vega(forward, strike, timeToExpiry, df, callOptionType)
        vegaPut = model.vega(forward, strike, timeToExpiry, df, putOptionType)

        assert round(vegaCall - vegaPut, dp) == 0.0, \
            "The method called 'vega()' doesn't comply with Call-Put parity"

        testCases.print("VEGA", vegaCall, vegaPut)

        #######################################################################

    except AssertionError as err:
        raise err

###############################################################################


test_Black()
testCases.compareTestCases()
