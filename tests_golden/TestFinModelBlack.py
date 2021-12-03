###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
# Guillaume Lefieux
###############################################################################

import numpy as np
from financepy.models.black import Black
from financepy.utils.global_types import OptionTypes
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)


def test_Black():

    forward = 0.034
    strike = 0.050
    riskFreeIR = 0.00
    time_to_expiry = 2.0
    volatility = 0.20

    testCases.header("ITEM", "CALL", "PUT")

    call_optionType = OptionTypes.EUROPEAN_CALL
    put_optionType = OptionTypes.EUROPEAN_PUT

    df = np.exp(-riskFreeIR * time_to_expiry)
    model = Black(volatility)

    dp = 12  # Precision

    try:

        #######################################################################

        valueCall = model.value(
            forward, strike, time_to_expiry, df, call_optionType)
        valuePut = model.value(
            forward, strike, time_to_expiry, df, put_optionType)

        assert round((valueCall - valuePut), dp) == round(df*(forward - strike), dp), \
            "The method called 'value()' doesn't comply with Call-Put parity"

        testCases.print("VALUE", valueCall, valuePut)

        #######################################################################

        deltaCall = model.delta(
            forward, strike, time_to_expiry, df, call_optionType)
        deltaPut = model.delta(
            forward, strike, time_to_expiry, df, put_optionType)

        assert round((1/df) * (deltaCall - deltaPut), dp) == 1.0, \
            "The method called 'delta()' doesn't comply with Call-put parity"

        testCases.print("DELTA", deltaCall, deltaPut)

        #######################################################################

        gammaCall = model.gamma(
            forward, strike, time_to_expiry, df, call_optionType)
        gammaPut = model.gamma(
            forward, strike, time_to_expiry, df, put_optionType)

        assert round(gammaCall - gammaPut, dp) == 0.0, \
            "The method called 'gamma()' doesn't comply with Call-Put parity"

        testCases.print("GAMMA", gammaCall, gammaPut)

        #######################################################################

        thetaCall = model.theta(
            forward, strike, time_to_expiry, df, call_optionType)
        thetaPut = model.theta(
            forward, strike, time_to_expiry, df, put_optionType)

        assert round((thetaCall - thetaPut), dp) == round((riskFreeIR * time_to_expiry) * (forward - strike) * df, dp), \
            "The method called 'theta()' doesn't comply with Call-Put parity"

        testCases.print("THETA", thetaCall, thetaPut)

        #######################################################################

        vegaCall = model.vega(
            forward, strike, time_to_expiry, df, call_optionType)
        vegaPut = model.vega(
            forward, strike, time_to_expiry, df, put_optionType)

        assert round(vegaCall - vegaPut, dp) == 0.0, \
            "The method called 'vega()' doesn't comply with Call-Put parity"

        testCases.print("VEGA", vegaCall, vegaPut)

        #######################################################################

    except AssertionError as err:
        raise err

###############################################################################


test_Black()
testCases.compareTestCases()
