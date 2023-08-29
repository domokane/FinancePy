###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys
import pandas as pd
sys.path.append("..")

from FinTestCases import FinTestCases, globalTestCaseMode
from financepy.products.bonds.bond_zero import BondZero
from financepy.products.bonds.bond import Bond, YTMCalcType

from financepy.products.bonds.zero_curve import BondZeroCurve

from financepy.utils.date import Date, from_datetime
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.math import ONE_MILLION

testCases = FinTestCases(__file__, globalTestCaseMode)

plotGraphs = False


###############################################################################

def test_bond_zero():
    issue_date = Date(25, 7, 2022)
    maturity_date = Date(24, 10, 2022)
    face_amount = 100.0
    issue_price = 99.6410

    bond = BondZero(issue_date=issue_date,
                    maturity_date=maturity_date,
                    issue_price=issue_price)

    settlement_date = Date(8, 8, 2022)

    clean_price = 99.6504

    ytm = bond.yield_to_maturity(settlement_date,
                                 clean_price,
                                 YTMCalcType.ZERO)

    accrued_interest = bond.accrued_interest(settlement_date, face_amount)

    testCases.header('YTM', 'accrued')
    testCases.print(ytm, accrued_interest)


def test_bond_zero_ror():
    test_case_file = 'test_cases_bond_zero_ror.csv'
    df = pd.read_csv('./data/' + test_case_file, parse_dates=['buy_date', 'sell_date'])
    # A 1-year bond with zero coupon per year. code: 092103011
    bond = BondZero(
        issue_date=Date(23, 7, 2021),
        maturity_date=Date(24, 8, 2022),
        issue_price=97.67
    )
    testCases.header('bond_code', 'buy_date', 'buy_ytm', 'buy_price', 'sell_date', 'sell_ytm', 'sell_price',
                     'simple_return', 'irr')
    for row in df.itertuples(index=False):
        buy_date = Date(row.buy_date.day, row.buy_date.month, row.buy_date.year)
        sell_date = Date(row.sell_date.day, row.sell_date.month, row.sell_date.year)
        buy_price = bond.dirty_price_from_ytm(buy_date, row.buy_ytm, YTMCalcType.ZERO)
        sell_price = bond.dirty_price_from_ytm(sell_date, row.sell_ytm, YTMCalcType.ZERO)
        simple, irr, pnl = bond.calc_ror(buy_date, sell_date, row.buy_ytm, row.sell_ytm)
        testCases.print(row.bond_code, buy_date, row.buy_ytm, buy_price, sell_date, row.sell_ytm, sell_price,
                        simple, irr)

###############################################################################


test_bond_zero()
test_bond_zero_ror()
testCases.compareTestCases()
