# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCount, DayCountTypes
from financepy.utils.date import Date


start = Date(1, 1, 2019)
end = Date(21, 5, 2019)
fin_freq = FrequencyTypes.ANNUAL

########################################################################################


def test_year_frace_thirty_360_bond():

    dc_type = DayCountTypes.THIRTY_360_BOND
    day_count = DayCount(dc_type)
    answer = day_count.year_frac(start, end, end, fin_freq)

    assert round(answer[0], 4) == 0.3889

########################################################################################


def test_year_frace_thirty_e_360():

    dc_type = DayCountTypes.THIRTY_E_360
    day_count = DayCount(dc_type)
    answer = day_count.year_frac(start, end, end, fin_freq)

    assert round(answer[0], 4) == 0.3889

########################################################################################


def test_year_frace_thirty_e_360_isda():

    dc_type = DayCountTypes.THIRTY_E_360_ISDA
    day_count = DayCount(dc_type)
    answer = day_count.year_frac(start, end, end, fin_freq)

    assert round(answer[0], 4) == 0.3889

########################################################################################


def test_year_frace_thirty_e_plus_360():

    dc_type = DayCountTypes.THIRTY_E_PLUS_360
    day_count = DayCount(dc_type)
    answer = day_count.year_frac(start, end, end, fin_freq)

    assert round(answer[0], 4) == 0.3889

########################################################################################


def test_year_frace_act_act_isda():

    dc_type = DayCountTypes.ACT_ACT_ISDA
    day_count = DayCount(dc_type)
    answer = day_count.year_frac(start, end, end, fin_freq)

    assert round(answer[0], 4) == 0.3836

########################################################################################


def test_year_frace_act_act_icma():

    dc_type = DayCountTypes.ACT_ACT_ICMA
    day_count = DayCount(dc_type)
    answer = day_count.year_frac(start, end, end, fin_freq)

    assert round(answer[0], 4) == 1.0000

########################################################################################


def test_year_frace_act_365_f():

    dc_type = DayCountTypes.ACT_365F
    day_count = DayCount(dc_type)
    answer = day_count.year_frac(start, end, end, fin_freq)

    assert round(answer[0], 4) == 0.3836

########################################################################################


def test_year_frace_act_360():

    dc_type = DayCountTypes.ACT_360
    day_count = DayCount(dc_type)
    answer = day_count.year_frac(start, end, end, fin_freq)

    assert round(answer[0], 4) == 0.3889

########################################################################################


def test_year_frace_act_365_l():

    dc_type = DayCountTypes.ACT_365L
    day_count = DayCount(dc_type)
    answer = day_count.year_frac(start, end, end, fin_freq)

    assert round(answer[0], 4) == 0.3836

########################################################################################


def test_year_frace_simple():

    dc_type = DayCountTypes.SIMPLE
    day_count = DayCount(dc_type)
    answer = day_count.year_frac(start, end, end, fin_freq)

    assert round(answer[0], 4) == 0.3836
