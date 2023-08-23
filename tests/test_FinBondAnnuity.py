###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.calendar import DateGenRuleTypes
from financepy.utils.calendar import BusDayAdjustTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.calendar import CalendarTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.date import Date, DateFormatTypes
from financepy.products.bonds.bond_annuity import BondAnnuity


def test_SemiAnnual_BondAnnuity():

    settlement_date = Date(20, 6, 2018)
    face = 1000000

    # Semi-Annual Frequency
    maturity_date = Date(20, 6, 2019)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
    basis_type = DayCountTypes.ACT_360

    annuity = BondAnnuity(maturity_date,
                          coupon,
                          freq_type,
                          calendar_type,
                          bus_day_adjust_type,
                          date_gen_rule_type,
                          basis_type)

    annuity.calculate_payments(settlement_date, face)

    assert len(annuity._flow_amounts) == 2 * 1 + 1
    assert len(annuity._coupon_dates) == 2 * 1 + 1

    assert annuity._coupon_dates[0] == settlement_date
    assert annuity._coupon_dates[-1] == maturity_date

    assert annuity._flow_amounts[0] == 0.0
    assert round(annuity._flow_amounts[-1]) == 25278.0

    assert annuity.accrued_interest(settlement_date, face) == 0.0


def test_Quarterly_BondAnnuity():

    settlement_date = Date(20, 6, 2018)
    face = 1000000

    # Quarterly Frequency
    maturity_date = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.QUARTERLY
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
    basis_type = DayCountTypes.ACT_360

    annuity = BondAnnuity(
        maturity_date,
        coupon,
        freq_type,
        calendar_type,
        bus_day_adjust_type,
        date_gen_rule_type,
        basis_type)

    annuity.calculate_payments(settlement_date, face)

    assert len(annuity._flow_amounts) == 10 * 4 + 1
    assert len(annuity._coupon_dates) == 10 * 4 + 1

    assert annuity._coupon_dates[0] == settlement_date
    assert annuity._coupon_dates[-1] == maturity_date

    assert annuity._flow_amounts[0] == 0.0
    assert round(annuity._flow_amounts[-1]) == 12778.0

    assert annuity.accrued_interest(settlement_date, face) == 0.0


def test_Monthly_BondAnnuity():

    settlement_date = Date(20, 6, 2018)
    face = 1000000

    # Monthly Frequency
    maturity_date = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.MONTHLY
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.BACKWARD
    basis_type = DayCountTypes.ACT_360

    annuity = BondAnnuity(maturity_date,
                          coupon,
                          freq_type,
                          calendar_type,
                          bus_day_adjust_type,
                          date_gen_rule_type,
                          basis_type)

    annuity.calculate_payments(settlement_date, face)

    assert len(annuity._flow_amounts) == 10*12 + 1
    assert len(annuity._coupon_dates) == 10*12 + 1

    assert annuity._coupon_dates[0] == settlement_date
    assert annuity._coupon_dates[-1] == maturity_date

    assert annuity._flow_amounts[0] == 0.0
    assert round(annuity._flow_amounts[-1]) == 4028.0

    assert annuity.accrued_interest(settlement_date, face) == 0.0


def test_ForwardGen_BondAnnuity():

    settlement_date = Date(20, 6, 2018)
    face = 1000000

    maturity_date = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.ANNUAL
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.FORWARD
    basis_type = DayCountTypes.ACT_360

    annuity = BondAnnuity(maturity_date,
                          coupon,
                          freq_type,
                          calendar_type,
                          bus_day_adjust_type,
                          date_gen_rule_type,
                          basis_type)

    annuity.calculate_payments(settlement_date, face)

    assert len(annuity._flow_amounts) == 10 * 1 + 1
    assert len(annuity._coupon_dates) == 10 * 1 + 1

    assert annuity._coupon_dates[0] == settlement_date
    assert annuity._coupon_dates[-1] == maturity_date

    assert round(annuity._flow_amounts[0]) == 0.0
    assert round(annuity._flow_amounts[-1]) == 50694.0

    assert annuity.accrued_interest(settlement_date, face) == 0.0


def test_ForwardGenWithLongEndStub_BondAnnuity():

    settlement_date = Date(20, 6, 2018)
    face = 1000000

    maturity_date = Date(20, 6, 2028)
    coupon = 0.05
    freq_type = FrequencyTypes.SEMI_ANNUAL
    calendar_type = CalendarTypes.WEEKEND
    bus_day_adjust_type = BusDayAdjustTypes.FOLLOWING
    date_gen_rule_type = DateGenRuleTypes.FORWARD
    basis_type = DayCountTypes.ACT_360

    annuity = BondAnnuity(maturity_date,
                          coupon,
                          freq_type,
                          calendar_type,
                          bus_day_adjust_type,
                          date_gen_rule_type,
                          basis_type)

    annuity.calculate_payments(settlement_date, face)

    assert len(annuity._flow_amounts) == 10 * 2 + 1
    assert len(annuity._coupon_dates) == 10 * 2 + 1

    assert annuity._coupon_dates[0] == settlement_date
    assert annuity._coupon_dates[-1] == maturity_date

    assert round(annuity._flow_amounts[0]) == 0.0
    assert round(annuity._flow_amounts[-1]) == 25417.0

    assert annuity.accrued_interest(settlement_date, face) == 0.0
