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

from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCount, DayCountTypes
from financepy.utils.date import Date



########################################################################################


def test_fin_day_count():

    print("DAY_COUNT_METHOD", "START", "END", "ALPHA")

    fin_freq = FrequencyTypes.ANNUAL

    for day_count_method in DayCountTypes:

        start_dt = Date(1, 1, 2019)
        next_dt = start_dt
        num_days = 20
        day_count = DayCount(day_count_method)

        for _ in range(0, num_days):
            next_dt = next_dt.add_days(7)
            dcf = day_count.year_frac(start_dt, next_dt, next_dt, fin_freq)

            print(str(day_count_method), str(start_dt), str(next_dt), dcf[0])


########################################################################################

test_fin_day_count()
