###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.date import Date
from financepy.products.bonds.bond_frn import BondFRN
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes


def test_bond_frn_1():
    # https://data.bloomberglp.com/bat/sites/3/2017/07/SF-2017_Paul-Fjeldsted.pdf
    issue_date = Date(10, 11, 2010)
    maturity_date = Date(10, 11, 2021)
    quoted_margin = 0.0025
    freq_type = FrequencyTypes.QUARTERLY
    accrual_type = DayCountTypes.THIRTY_E_360

    bond = BondFRN(issue_date,
                   maturity_date,
                   quoted_margin,
                   freq_type,
                   accrual_type)

    clean_price = 96.793
    resetIbor = 0.0143456 - quoted_margin
    current_ibor = 0.0120534
    future_ibors = 0.0130522

    settlement_date = Date(21, 7, 2017)

    dm = bond.discount_margin(settlement_date,
                              resetIbor,
                              current_ibor,
                              future_ibors,
                              clean_price)

    assert round(dm * 10000, 4) == 103.1985

    dirty_price = bond.dirty_price_from_dm(settlement_date,
                                           resetIbor,
                                           current_ibor,
                                           future_ibors,
                                           dm)

    assert round(dirty_price, 4) == 97.0266

    lastCouponDt = bond._pcd
    assert lastCouponDt == Date(10, 5, 2017)

    accddays = bond._accrued_days
    assert accddays == 71

    accdAmount = bond._accrued_interest
    assert round(accdAmount, 4) == 0.0023

    principal = bond.principal(settlement_date,
                               resetIbor,
                               current_ibor,
                               future_ibors,
                               dm)

    assert round(principal, 4) == 97.0243

    duration = bond.dollar_duration(settlement_date,
                                    resetIbor,
                                    current_ibor,
                                    future_ibors,
                                    dm)

    assert round(duration, 4) == 5.1148

    modified_duration = bond.modified_duration(settlement_date,
                                               resetIbor,
                                               current_ibor,
                                               future_ibors,
                                               dm)

    assert round(modified_duration, 4) == 0.0527

    macauley_duration = bond.macauley_duration(settlement_date,
                                               resetIbor,
                                               current_ibor,
                                               future_ibors,
                                               dm)

    assert round(macauley_duration, 4) == 0.0530

    convexity = bond.convexity_from_dm(settlement_date,
                                       resetIbor,
                                       current_ibor,
                                       future_ibors,
                                       dm)

    assert round(convexity, 8) == 0.00005558

    duration = bond.dollar_credit_duration(settlement_date,
                                           resetIbor,
                                           current_ibor,
                                           future_ibors,
                                           dm)

    assert round(duration, 4) == 401.0636

    modified_duration = bond.modified_credit_duration(settlement_date,
                                                      resetIbor,
                                                      current_ibor,
                                                      future_ibors,
                                                      dm)

    assert round(modified_duration, 4) == 4.1335


def test_bond_frn_2():
    # https://ebrary.net/14293/economics/actual_floater
    issue_date = Date(28, 3, 2000)
    settlement_date = Date(28, 3, 2014)
    maturity_date = Date(3, 2, 2021)
    quoted_margin = 0.0020
    freq_type = FrequencyTypes.SEMI_ANNUAL
    accrual_type = DayCountTypes.THIRTY_E_360_ISDA

    bond = BondFRN(issue_date,
                   maturity_date,
                   quoted_margin,
                   freq_type,
                   accrual_type)

    clean_price = 93.08
    resetIbor = 0.00537 - quoted_margin
    current_ibor = 0.027558
    future_ibors = 0.03295

    dm = bond.discount_margin(settlement_date,
                              resetIbor,
                              current_ibor,
                              future_ibors,
                              clean_price)

    assert round(dm * 10000, 4) == 123.0623

    dirty_price = bond.dirty_price_from_dm(settlement_date,
                                           resetIbor,
                                           current_ibor,
                                           future_ibors,
                                           dm)

    assert round(dirty_price, 4) == 93.1315

    lastCouponDt = bond._pcd
    assert lastCouponDt == Date(3, 2, 2014)

    accddays = bond._accrued_days
    assert accddays == 55

    accdAmount = bond._accrued_interest
    assert round(accdAmount, 4) == 0.0005

    principal = bond.principal(settlement_date,
                               resetIbor,
                               current_ibor,
                               future_ibors,
                               dm)

    assert round(principal, 4) == 93.131

    duration = bond.dollar_duration(settlement_date,
                                    resetIbor,
                                    current_ibor,
                                    future_ibors,
                                    dm)

    assert round(duration, 4) == 31.8958

    modified_duration = bond.modified_duration(settlement_date,
                                               resetIbor,
                                               current_ibor,
                                               future_ibors,
                                               dm)

    assert round(modified_duration, 4) == 0.3425

    macauley_duration = bond.macauley_duration(settlement_date,
                                               resetIbor,
                                               current_ibor,
                                               future_ibors,
                                               dm)

    assert round(macauley_duration, 4) == 0.3452

    convexity = bond.convexity_from_dm(settlement_date,
                                       resetIbor,
                                       current_ibor,
                                       future_ibors,
                                       dm)

    assert round(convexity, 4) == 0.0023

    principal = bond.principal(settlement_date,
                               resetIbor,
                               current_ibor,
                               future_ibors,
                               dm)

    assert round(principal, 4) == 93.1310

    duration = bond.dollar_credit_duration(settlement_date,
                                           resetIbor,
                                           current_ibor,
                                           future_ibors,
                                           dm)

    assert round(duration, 4) == 563.2624

    modified_duration = bond.modified_credit_duration(settlement_date,
                                                      resetIbor,
                                                      current_ibor,
                                                      future_ibors,
                                                      dm)

    assert round(modified_duration, 4) == 6.0480
