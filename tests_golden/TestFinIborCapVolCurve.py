###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date
from financepy.market.volatility.ibor_cap_vol_curve import IborCapVolCurve
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

##########################################################################


def test_FinCapVolCurve():

    # Reproduces example in Table 32.1 of Hull Book
    valuation_date = Date(1, 1, 2020)

    capVolDates = []
    capletVolTenor = "1Y"
    num_periods = 10
    capletDt = valuation_date

    capVolDates.append(valuation_date)

    for i in range(0, num_periods):
        capletDt = capletDt.add_tenor(capletVolTenor)
        capVolDates.append(capletDt)

    capVolatilities = [0.0, 15.50, 18.25, 17.91, 17.74, 17.27,
                       16.79, 16.30, 16.01, 15.76, 15.54]
    capVolatilities = np.array(capVolatilities)/100.0

    day_count_type = DayCountTypes.ACT_ACT_ISDA
    volCurve = IborCapVolCurve(valuation_date,
                               capVolDates,
                               capVolatilities,
                               day_count_type)

    testCases.header("DATE", "CAPVOL", "CAPLETVOL")
    for dt in capVolDates:
        capFloorVol = volCurve.cap_vol(dt)
        capFloorLetVol = volCurve.caplet_vol(dt)
        testCases.print("%s" % dt,
                        "%7.3f" % (capFloorVol*100.0),
                        "%7.2f" % (capFloorLetVol*100.0))

##########################################################################


test_FinCapVolCurve()
testCases.compareTestCases()
