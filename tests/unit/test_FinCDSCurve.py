# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from financepy.products.credit.cds import CDS
from financepy.products.rates.ibor_swap import IborSwap
from financepy.market.curves.cds_curve import CDSCurve
from financepy.market.curves.ibor_single_curve import IborSingleCurve
from financepy.utils.frequency import FrequencyTypes
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date
from financepy.utils.global_types import SwapTypes
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


########################################################################################


import numpy as np


def test_fin_cds_curve():

    curve_dt = Date(20, 12, 2018)

    swaps = []
    depos = []
    fras = []

    fixed_dcc = DayCountTypes.ACT_365F
    fixed_freq = FrequencyTypes.SEMI_ANNUAL
    fixed_cpn = 0.05

    for i in range(1, 11):

        maturity_dt = curve_dt.add_months(12 * i)
        swap = IborSwap(
            curve_dt,
            maturity_dt,
            SwapTypes.PAY,
            fixed_cpn,
            fixed_freq,
            fixed_dcc,
        )
        swaps.append(swap)

    libor_curve = IborSingleCurve(curve_dt, depos, fras, swaps)

    cds_contracts = []

    for i in range(1, 11):
        maturity_dt = curve_dt.add_months(12 * i)
        cds = CDS(curve_dt, maturity_dt, 0.005 + 0.001 * (i - 1))
        cds_contracts.append(cds)

    recovery_rate = 0.40

    issuer_curve = CDSCurve(
        curve_dt, cds_contracts, libor_curve, recovery_rate
    )

    assert round(issuer_curve._times[0], 4) == 0.0
    assert round(issuer_curve._times[5], 4) == 5.0027
    assert round(issuer_curve._times[9], 4) == 9.0055
    assert round(issuer_curve._qs[0], 4) == 1.0
    assert round(issuer_curve._qs[5], 4) == 0.9249
    assert round(issuer_curve._qs[9], 4) == 0.8071

    for i in [1, 5, 10]:
        maturity_dt = curve_dt.add_months(12 * i)
        cds = CDS(curve_dt, maturity_dt, 0.005 + 0.001 * (i - 1))
        dirty_pv, clean_pv = cds.value(curve_dt, issuer_curve, recovery_rate)
        assert abs(dirty_pv) < 1e-5
        assert abs(clean_pv) < 1e-5

    # Test scalar Date evaluation and day-count consistency
    d1 = curve_dt.add_months(6)
    d2 = curve_dt.add_months(12)
    df_single = issuer_curve.df(d1)
    sp_single = issuer_curve.survival_prob(d1)

    assert isinstance(df_single, float)
    assert isinstance(sp_single, float)
    assert df_single == libor_curve.df(d1)

    # Test vectorization with list of Dates
    dts = [d1, d2]
    df_vector = issuer_curve.df(dts)
    sp_vector = issuer_curve.survival_prob(dts)

    assert isinstance(df_vector, np.ndarray)
    assert isinstance(sp_vector, np.ndarray)
    assert len(df_vector) == 2
    assert len(sp_vector) == 2
    assert np.isclose(df_vector[0], df_single)
    assert np.isclose(df_vector[1], issuer_curve.df(d2))
    assert np.isclose(sp_vector[0], sp_single)
    assert np.isclose(sp_vector[1], issuer_curve.survival_prob(d2))

    # Test numeric float and numpy array inputs
    t_single = 0.5
    df_t = issuer_curve.df(t_single)
    sp_t = issuer_curve.survival_prob(t_single)
    assert isinstance(df_t, float)
    assert isinstance(sp_t, float)
    assert np.isclose(df_t, libor_curve.df_t(t_single))

    t_arr = np.array([0.5, 1.0])
    df_arr = issuer_curve.df(t_arr)
    sp_arr = issuer_curve.survival_prob(t_arr)
    assert isinstance(df_arr, np.ndarray)
    assert isinstance(sp_arr, np.ndarray)
    assert len(df_arr) == 2
    assert len(sp_arr) == 2
    assert np.isclose(df_arr[0], df_t)
    assert np.isclose(sp_arr[0], sp_t)

    # Test day count consistency when libor_curve uses a different day count
    libor_curve_360 = IborSingleCurve(
        curve_dt, depos, fras, swaps, curve_dc_type=DayCountTypes.ACT_360
    )
    issuer_curve_360 = CDSCurve(
        curve_dt, cds_contracts, libor_curve_360, recovery_rate
    )
    assert issuer_curve_360.df(d1) == libor_curve_360.df(d1)
    assert issuer_curve_360.df(dts)[0] == libor_curve_360.df(dts)[0]
    assert issuer_curve_360.df(dts)[1] == libor_curve_360.df(dts)[1]
