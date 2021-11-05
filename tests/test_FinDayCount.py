###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCount, DayCountTypes
from financepy.utils.date import Date


start = Date(1, 1, 2019)
end = Date(21, 5, 2019)
finFreq = FrequencyTypes.ANNUAL


def test_year_frace_THIRTY_360_BOND():
    day_count_type = DayCountTypes.THIRTY_360_BOND
    day_count = DayCount(day_count_type)
    answer = day_count.year_frac(start, end, end, finFreq)

    assert round(answer[0], 4) == 0.3889


def test_year_frace_THIRTY_E_360():
    day_count_type = DayCountTypes.THIRTY_E_360
    day_count = DayCount(day_count_type)
    answer = day_count.year_frac(start, end, end, finFreq)

    assert round(answer[0], 4) == 0.3889


def test_year_frace_THIRTY_E_360_ISDA():
    day_count_type = DayCountTypes.THIRTY_E_360_ISDA
    day_count = DayCount(day_count_type)
    answer = day_count.year_frac(start, end, end, finFreq)

    assert round(answer[0], 4) == 0.3889


def test_year_frace_THIRTY_E_PLUS_360():
    day_count_type = DayCountTypes.THIRTY_E_PLUS_360
    day_count = DayCount(day_count_type)
    answer = day_count.year_frac(start, end, end, finFreq)

    assert round(answer[0], 4) == 0.3889


def test_year_frace_ACT_ACT_ISDA():
    day_count_type = DayCountTypes.ACT_ACT_ISDA
    day_count = DayCount(day_count_type)
    answer = day_count.year_frac(start, end, end, finFreq)

    assert round(answer[0], 4) == 0.3836


def test_year_frace_ACT_ACT_ICMA():
    day_count_type = DayCountTypes.ACT_ACT_ICMA
    day_count = DayCount(day_count_type)
    answer = day_count.year_frac(start, end, end, finFreq)

    assert round(answer[0], 4) == 1.0000


def test_year_frace_ACT_365F():
    day_count_type = DayCountTypes.ACT_365F
    day_count = DayCount(day_count_type)
    answer = day_count.year_frac(start, end, end, finFreq)

    assert round(answer[0], 4) == 0.3836


def test_year_frace_ACT_360():
    day_count_type = DayCountTypes.ACT_360
    day_count = DayCount(day_count_type)
    answer = day_count.year_frac(start, end, end, finFreq)

    assert round(answer[0], 4) == 0.3889


def test_year_frace_ACT_365L():
    day_count_type = DayCountTypes.ACT_365L
    day_count = DayCount(day_count_type)
    answer = day_count.year_frac(start, end, end, finFreq)

    assert round(answer[0], 4) == 0.3836


def test_year_frace_SIMPLE():
    day_count_type = DayCountTypes.SIMPLE
    day_count = DayCount(day_count_type)
    answer = day_count.year_frac(start, end, end, finFreq)

    assert round(answer[0], 4) == 0.3836
