# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.
import sys as _sys
from pathlib import Path as _Path

_EXAMPLES_CODE = _Path(__file__).resolve().parents[1]
if str(_EXAMPLES_CODE) not in _sys.path:
    _sys.path.insert(0, str(_EXAMPLES_CODE))
from double_click_pause import install_double_click_pause as _install_double_click_pause

_install_double_click_pause()
import add_fp_to_path

from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.schedule import Schedule
from financepy.utils.date import Date

########################################################################################


def test_dt_adjust():

    start_dt = Date(28, 2, 2008)
    end_dt = Date(28, 2, 2011)

    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.NONE
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    print("NO ADJUSTMENTS", "DATE")
    schedule = Schedule(start_dt, end_dt, freq_type, cal_type, bd_type, dg_type)

    for dt in schedule.adjusted_dts:
        print("Date:", dt)

    print("")
    print("NO WEEKENDS AND FOLLOWING", "DATE")
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    schedule = Schedule(start_dt, end_dt, freq_type, cal_type, bd_type, dg_type)

    for dt in schedule.adjusted_dts:
        print("Date:", dt)

    print("")
    print("NO WEEKENDS AND MODIFIED FOLLOWING", "DATE")
    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.WEEKEND
    bd_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    schedule = Schedule(start_dt, end_dt, freq_type, cal_type, bd_type, dg_type)

    for dt in schedule.adjusted_dts:
        print("Date:", dt)

    print("")
    print("NO WEEKENDS AND US HOLIDAYS AND MODIFIED FOLLOWING", "DATE")

    freq_type = FrequencyTypes.SEMI_ANNUAL
    cal_type = CalendarTypes.UNITED_STATES
    bd_type = BusDayAdjustTypes.MODIFIED_FOLLOWING
    dg_type = DateGenRuleTypes.BACKWARD

    start_dt = Date(4, 7, 2008)
    end_dt = Date(4, 7, 2011)

    schedule = Schedule(start_dt, end_dt, freq_type, cal_type, bd_type, dg_type)

    for dt in schedule.adjusted_dts:
        print("Date:", dt)


########################################################################################

test_dt_adjust()
