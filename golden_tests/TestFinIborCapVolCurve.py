###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np

import sys

sys.path.append("..")

from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date
from financepy.market.volatility.ibor_cap_vol_curve import IborCapVolCurve
from FinTestCases import FinTestCases, globalTestCaseMode

test_cases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_FinCapVolCurve():

    # Reproduces example in Table 32.1 of Hull Book
    value_dt = Date(1, 1, 2020)

    capVolDates = []
    capletVolTenor = "1Y"
    num_periods = 10
    capletDt = value_dt

    capVolDates.append(value_dt)

    for i in range(0, num_periods):
        capletDt = capletDt.add_tenor(capletVolTenor)
        capVolDates.append(capletDt)

    capVolatilities = [
        0.0,
        15.50,
        18.25,
        17.91,
        17.74,
        17.27,
        16.79,
        16.30,
        16.01,
        15.76,
        15.54,
    ]
    capVolatilities = np.array(capVolatilities) / 100.0

    day_count_type = DayCountTypes.ACT_ACT_ISDA
    vol_curve = IborCapVolCurve(
        value_dt, capVolDates, capVolatilities, day_count_type
    )

    test_cases.header("DATE", "CAPVOL", "CAPLETVOL")
    for dt in capVolDates:
        capFloorVol = vol_curve.cap_vol(dt)
        capFloorLetVol = vol_curve.caplet_vol(dt)
        test_cases.print(
            "%s" % dt,
            "%7.3f" % (capFloorVol * 100.0),
            "%7.2f" % (capFloorLetVol * 100.0),
        )


##########################################################################


test_FinCapVolCurve()
test_cases.compareTestCases()
