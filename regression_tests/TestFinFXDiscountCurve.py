# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import add_fp_to_path

from financepy.utils.date import Date
from financepy.utils.frequency import FrequencyTypes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.market.curves.fx_discount_curve import FXDiscountCurve
from FinTestCases import FinTestCases, global_test_case_mode

test_cases = FinTestCases(__file__, global_test_case_mode)


########################################################################################
def test_fin_fx_discount_curve():
    value_dt = Date(1, 1, 2024)

    spot_fx_rate = 150.0
    forward_tenors = ["1M", "3M", "6M"]
    forward_dts = value_dt.add_months([1, 3, 6])
    forward_fx_rates = [149.8, 149.2, 148.5]

    domestic_curve = DiscountCurveFlat(
        value_dt,
        0.05,
        FrequencyTypes.CONTINUOUS,
    )

    foreign_curve = FXDiscountCurve(
        value_dt,
        spot_fx_rate,
        forward_dts,
        forward_fx_rates,
        domestic_curve,
    )

    test_cases.header("TENOR", "INPUT_FWD", "IMPLIED_FWD")

    for tenor, forward_dt, forward_fx_rate in zip(
        forward_tenors,
        forward_dts,
        forward_fx_rates,
    ):
        implied_forward = (
            spot_fx_rate
            * foreign_curve.df(forward_dt)
            / domestic_curve.df(forward_dt)
        )

        test_cases.print(tenor, forward_fx_rate, implied_forward)


########################################################################################

test_fin_fx_discount_curve()
test_cases.compare_test_cases()
