###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import sys

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

def testBondZeroCoupon():
    issue_date = Date(25, 7, 2022)
    maturity_date = Date(24, 10, 2022)
    face_amount = ONE_MILLION
    issue_price = 99.6410

    bond = BondZero(issue_date=issue_date,
                    maturity_date=maturity_date,
                    face_amount=face_amount,
                    issue_price=issue_price)

    settlement_date = Date(8, 8, 2022)

    clean_price = 99.6504
    ytm = bond.yield_to_maturity(settlement_date,
                                 clean_price,
                                 YTMCalcType.ZERO)
    accrued_interest = bond.calc_accrued_interest(settlement_date)

    testCases.header('YTM', 'accrued')
    testCases.print(ytm, accrued_interest)


###############################################################################

testBondZeroCoupon()
testCases.compareTestCases()
