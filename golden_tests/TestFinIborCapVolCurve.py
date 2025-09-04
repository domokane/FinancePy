# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import numpy as np

import add_fp_to_path

from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date
from financepy.market.volatility.ibor_cap_vol_curve import IborCapVolCurve
from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)

########################################################################################


def test_fin_cap_vol_curve():

    # Reproduces example in Table 32.1 of Hull Book
    value_dt = Date(1, 1, 2020)

    cap_vol_dates = []
    caplet_vol_tenor = "1Y"
    num_periods = 10
    caplet_dt = value_dt

    cap_vol_dates.append(value_dt)

    for i in range(0, num_periods):
        caplet_dt = caplet_dt.add_tenor(caplet_vol_tenor)
        cap_vol_dates.append(caplet_dt)

    cap_volatilities = [
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
    cap_volatilities = np.array(cap_volatilities) / 100.0

    dc_type = DayCountTypes.ACT_ACT_ISDA
    vol_curve = IborCapVolCurve(value_dt, cap_vol_dates, cap_volatilities, dc_type)

    test_cases.header("DATE", "CAPVOL", "CAPLETVOL")
    for dt in cap_vol_dates:
        cap_floor_vol = vol_curve.cap_vol(dt)
        cap_floor_let_vol = vol_curve.caplet_vol(dt)
        test_cases.print(
            "%s" % dt,
            "%7.3f" % (cap_floor_vol * 100.0),
            "%7.2f" % (cap_floor_let_vol * 100.0),
        )


########################################################################################

test_fin_cap_vol_curve()
test_cases.compare_test_cases()
