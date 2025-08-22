########################################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
########################################################################################

import sys, os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from .helpers import build_ibor_curve, build_issuer_curve
from financepy.utils.global_types import SwapTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.credit.cds_curve import CDSCurve
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.rates.ibor_swap import IborSwap
from financepy.utils.math import ONE_MILLION
from financepy.products.credit.cds import CDS


# We treat an index as a CDS contract with a flat CDS curve
trade_dt = Date(7, 2, 2006)
libor_curve = build_ibor_curve(trade_dt)
issuer_curve = build_issuer_curve(trade_dt, libor_curve)
step_in_dt = trade_dt.add_days(1)
value_dt = step_in_dt
maturity_dt = Date(20, 6, 2010)

cdsRecovery = 0.40
notional = 10.0 * ONE_MILLION
long_protection = True
index_cpn = 0.004

cds_indexContract = CDS(step_in_dt, maturity_dt, index_cpn, notional, long_protection)


def test_cds_index():
    spd = cds_indexContract.par_spread(value_dt, issuer_curve, cdsRecovery) * 10000.0
    assert round(spd, 4) == 48.3748

    v = cds_indexContract.value(value_dt, issuer_curve, cdsRecovery)
    assert round(v["dirty_pv"], 4) == 27064.9906
    assert round(v["clean_pv"], 4) == 32620.5461

    p = cds_indexContract.clean_price(value_dt, issuer_curve, cdsRecovery)
    assert round(p, 4) == 99.6738

    accrued_days = cds_indexContract.accrued_days()
    assert accrued_days == 50.0

    accrued_interest = cds_indexContract.accrued_interest()
    assert round(accrued_interest, 4) == -5555.5556

    prot_pv = cds_indexContract.prot_leg_pv(value_dt, issuer_curve, cdsRecovery)
    assert round(prot_pv, 4) == 188423.9948

    premPV = cds_indexContract.premium_leg_pv(value_dt, issuer_curve, cdsRecovery)
    assert round(premPV, 4) == 161359.0042
