# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.utils.global_types import SwapTypes
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.products.credit.cds_curve import CDSCurve
from financepy.products.rates.ibor_single_curve import IborSingleCurve
from financepy.products.rates.ibor_swap import IborSwap
from financepy.utils.math import ONE_MILLION
from financepy.products.credit.cds import CDS


from .helpers import build_ibor_curve, build_issuer_curve

# We treat an index as a CDS contract with a flat CDS curve
trade_dt = Date(7, 2, 2006)
libor_curve = build_ibor_curve(trade_dt)
issuer_curve = build_issuer_curve(trade_dt, libor_curve)
step_in_dt = trade_dt.add_days(1)
value_dt = step_in_dt
maturity_dt = Date(20, 6, 2010)

cds_recovery = 0.40
notional = 10.0 * ONE_MILLION
long_protection = True
index_cpn = 0.004

cds_index_contract = CDS(step_in_dt, maturity_dt, index_cpn, notional, long_protection)

########################################################################################


def test_cds_index():

    spd = cds_index_contract.par_spread(value_dt, issuer_curve, cds_recovery) * 10000.0
    assert round(spd, 4) == 48.3748

    v = cds_index_contract.value(value_dt, issuer_curve, cds_recovery)
    assert round(v["dirty_pv"], 4) == 27063.9928
    assert round(v["clean_pv"], 4) == 32619.5484

    p = cds_index_contract.clean_price(value_dt, issuer_curve, cds_recovery)
    assert round(p, 4) == 99.6738

    accrued_days = cds_index_contract.accrued_days()
    assert accrued_days == 50.0

    accrued_interest = cds_index_contract.accrued_interest()
    assert round(accrued_interest, 4) == -5555.5556

    prot_pv = cds_index_contract.prot_leg_pv(value_dt, issuer_curve, cds_recovery)
    assert round(prot_pv, 4) == 188418.0071

    prem_pv = cds_index_contract.premium_leg_pv(value_dt, issuer_curve, cds_recovery)
    assert round(prem_pv, 4) == 161354.0143
