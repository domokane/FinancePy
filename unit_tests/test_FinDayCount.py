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


def test_act_365_l_without_reference_end_date():

    day_count = DayCount(DayCountTypes.ACT_365L)
    period_start = Date(1, 12, 2023)
    period_end = Date(1, 3, 2024)

    implicit_period_end = day_count.year_frac(
        period_start,
        period_end,
        freq_type=FrequencyTypes.ANNUAL,
    )
    explicit_period_end = day_count.year_frac(
        period_start,
        period_end,
        period_end,
        FrequencyTypes.ANNUAL,
    )

    assert implicit_period_end == explicit_period_end
    assert implicit_period_end == (91.0 / 366.0, 91.0, 366)

    # A separately supplied coupon-period end must still control the
    # denominator for an accrued fraction.
    accrued_start = Date(1, 3, 2023)
    accrued_end = Date(1, 12, 2023)
    next_coupon = Date(1, 3, 2024)

    full_period_mode = day_count.year_frac(
        accrued_start,
        accrued_end,
        freq_type=FrequencyTypes.ANNUAL,
    )
    accrued_mode = day_count.year_frac(
        accrued_start,
        accrued_end,
        next_coupon,
        FrequencyTypes.ANNUAL,
    )

    assert full_period_mode == (275.0 / 365.0, 275.0, 365)
    assert accrued_mode == (275.0 / 366.0, 275.0, 366)


########################################################################################


def test_year_frace_simple():

    dc_type = DayCountTypes.JULIAN
    day_count = DayCount(dc_type)
    answer = day_count.year_frac(start, end, end, fin_freq)

    # 1 Jan 2019 -> 21 May 2019 is 140 days; 140 / 365.25 = 0.38330
    assert round(answer[0], 4) == 0.3833

    # A Julian year is 365.25 days, so this must not collapse onto ACT_365F.
    act_365f = DayCount(DayCountTypes.ACT_365F).year_frac(start, end, end, fin_freq)
    assert answer[0] != act_365f[0]
