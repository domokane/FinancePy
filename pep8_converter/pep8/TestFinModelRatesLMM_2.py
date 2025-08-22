# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from FinTestCases import FinTestCases, global_test_case_mode
from financepy.models.lmm_mc import lmm_sticky_caplet_pricer
from financepy.models.lmm_mc import lmm_ratchet_caplet_pricer
from financepy.models.lmm_mc import lmm_fwd_fwd_correlation
from financepy.models.lmm_mc import lmm_swap_pricer
from financepy.models.lmm_mc import lmm_price_caps_black
from financepy.models.lmm_mc import lmm_cap_flr_pricer
from financepy.models.lmm_mc import lmm_swaption_vol_approx
from financepy.models.lmm_mc import lmm_sim_swaption_vol
from financepy.models.lmm_mc import lmm_swaption_pricer
from financepy.models.lmm_mc import lmm_simulate_fwds_mf
from financepy.models.lmm_mc import lmm_simulate_fwds_1f
from financepy.models.lmm_mc import lmm_simulate_fwds_nf
from financepy.utils.helpers import check_vector_differences
from financepy.products.rates.ibor_swaption import IborSwaption
from financepy.products.rates.ibor_swaption import SwapTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black import Black
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date
from financepy.market.volatility.ibor_cap_vol_curve import IborCapVolCurve
import numpy as np
import time as time
import sys

sys.path.append("..")


test_cases = FinTestCases(__file__, global_test_case_mode)



########################################################################################


def get_correlation_matrix(num_fwds, beta, dt):
    correl = np.zeros((num_fwds, num_fwds))  # indexes from 1 to p-1
    for i in range(0, num_fwds):
        correl[i][i] = 1.0
        for j in range(i + 1, num_fwds):
            correl[i][j] = np.exp(-beta * (j - i) * dt)
            correl[j][i] = correl[i][j]

    return correl


########################################################################################


""" def get_vol_curve(num_fwds, dt, flat_vol=None):

    value_dt = Date(1, 1, 2020)

    cap_vol_dates = []
    caplet_vol_tenor = "1Y"
    num_periods = 10
    caplet_dt = value_dt

    cap_vol_dates.append(value_dt)
    for _ in range(0, num_periods):
        caplet_dt = caplet_dt.add_tenor(caplet_vol_tenor)
        cap_vol_dates.append(caplet_dt)

    if flat_vol is None:
        cap_volatilities = [0.0, 15.50, 18.25, 17.91, 17.74, 17.27,
                           16.79, 16.30, 16.01, 15.76, 15.54]
        cap_volatilities = np.array(cap_volatilities)/100.0
    else:
        cap_volatilities = [flat_vol] * (num_periods+1)
        cap_volatilities = np.array(cap_volatilities)
        cap_volatilities[0] = 0.0

    dc_type = DayCountTypes.ACT_ACT_ISDA
    vol_curve = IborCapVolCurve(value_dt,
                                   cap_vol_dates,
                                   cap_volatilities,
                                   dc_type)

    zetas = np.zeros(num_fwds)
    t = 0.0
    for ix in range(0, num_fwds):
        t = t + dt
        zetas[ix] = vol_curve.cap_floorlet_vol(t)

    print(zetas)
    return zetas """

########################################################################################


########################################################################################


def get_forward_curve(num_fwds, r):

    fwd0 = np.zeros(num_fwds)
    for i in range(0, num_fwds):
        fwd0[i] = r

    fwd0 = np.array(fwd0)
    return fwd0


########################################################################################


# def test_Swaptions():

#     dt = 0.5
#     t_exp = 3.0
#     t_mat = 10.0
#     a = int(2*t_exp)
#     b = int(2*t_mat)
#     num_fwds = 20
#     taus = np.array([dt] * num_fwds)

#     r = 0.05
#     fwd0 = get_forward_curve(num_fwds, r)
#     correl = get_correlation_matrix(num_fwds, 0.00000000000001, dt)

#     fwd_rateVol = 0.20
#     zetas = get_vol_curve(num_fwds, dt, fwd_rateVol)

#     seed = 1489
#     num_paths = 2000  # 100000
#     fwds_nf = lmm_simulate_fwds_nf(num_fwds, num_paths, fwd0,
#                                zetas, correl, taus, seed)
#     strike = r
#     PAYSwaption = 1
#     use_sobol = 0
#     numeraire_index = 0

#     fwds1F = lmm_simulate_fwds1_f(num_fwds, num_paths, numeraire_index, fwd0,
#                                zetas, taus, use_sobol, seed)

#     for i_exp in range(1, 10):

#         t_exp = float(i_exp)
#         a = int(2*t_exp)
#         print(a, b)

#         swaption_price1F = lmm_swaption_pricer(strike, a, b, num_paths,
#                                             fwd0, fwds1F, taus, PAYSwaption)

#         swaption_priceNF = lmm_swaption_pricer(strike, a, b, num_paths,
#                                             fwd0, fwds_nf, taus, PAYSwaption)

#         swaption_vol = lmm_swaption_vol_approx(a, b, fwd0, taus, zetas, correl)

#         swap_vol_sim1_f = lmm_sim_swaption_vol(a, b, fwd0, fwds1F, taus)
#         swap_vol_sim_nf = lmm_sim_swaption_vol(a, b, fwd0, fwds_nf, taus)

#         value_dt = Date(1, 1, 2010)
#         libor_curve = fin_discount_curve_flat(value_dt, r,
#                                           FrequencyTypes.QUARTERLY)

#         settle_dt = value_dt
#         exercise_dt = settle_dt.add_months(a*3)
#         maturity_dt = settle_dt.add_months(b*3)

#         fixed_cpn = strike
#         fixed_freq_type = FrequencyTypes.QUARTERLY
#         fixed_dc_type = DayCountTypes.ACT_ACT_ISDA
#         float_freq_type = FrequencyTypes.QUARTERLY
#         float_dc_type = DayCountTypes.ACT_ACT_ISDA
#         notional = 1.0

#         # Pricing a PAY
#         swaption_type = ibor_swaption_types.PAY
#         swaption = IborSwaption(settle_dt,
#                                     exercise_dt,
#                                     maturity_dt,
#                                     swaption_type,
#                                     fixed_cpn,
#                                     fixed_freq_type,
#                                     fixed_dc_type,
#                                     notional,
#                                     float_freq_type,
#                                     float_dc_type)

#         model = Black(swaption_vol)
#         black_swaption_price = swaption.value(value_dt, libor_curve, model)

#         print("K:%6.5f t_exp:%8.2f fwd_vol:%9.5f sim_vol1_f:%9.5f sim_vol_nf:%9.5f reb_vol:%9.5f sim_px1_f:%9.5f sim_px_nf:%9.5f Black Px:%9.5f"
#               % (strike, t_exp, fwd_rateVol, swap_vol_sim1_f, swap_vol_sim_nf, swaption_vol,
#                  swaption_price1F, swaption_priceNF, black_swaption_price))

# #        print(swaption)

########################################################################################


# def test_CapsFloors():

#     num_fwds = 40
#     num_paths = 100000
#     dt = 0.25
#     taus = np.array([dt] * num_fwds)
#     seeds = [42] #, 143, 101, 993, 9988]

#     r = 0.04
#     fwd0 = get_forward_curve(num_fwds, r)

#     fwd_rateVol = 0.20
#     zetas = get_vol_curve(num_fwds, dt, fwd_rateVol)

#     correl = get_correlation_matrix(num_fwds, 100.0, dt)

#     # At the money
#     K = r
#     caplet_prices_black = lmm_price_caps_black(fwd0, zetas, num_fwds, K, taus)

#     num_factors = 1
#     numeraire_index = 1
#     use_sobol = 1

#     # Examine variance for different seeds
#     for seed in seeds:

#         print("=============================================================")
#         print("Seed:", seed)

#         fwds1F = lmm_simulate_fwds1_f(num_fwds, num_paths, numeraire_index, fwd0,
#                                    zetas, taus, use_sobol, seed)

#         sum_cap1_f = lmm_cap_flr_pricer(num_fwds, num_paths, K, fwd0, fwds1F, taus, 1)
#         sum_flr1_f = lmm_cap_flr_pricer(num_fwds, num_paths, K, fwd0, fwds1F, taus, 0)

#         fwds_nf = lmm_simulate_fwds_nf(num_fwds, num_paths, fwd0,
#                                    zetas, correl, taus, seed)

#         sum_cap_nf = lmm_cap_flr_pricer(num_fwds, num_paths, K, fwd0, fwds_nf, taus, 1)
#         sum_flr_nf = lmm_cap_flr_pricer(num_fwds, num_paths, K, fwd0, fwds_nf, taus, 0)

#         print("i     CAPLET1F    FLRLET1F   CAPLETNF    FLRLET_NF   BS CAP")
#         for i in range(0, num_fwds):
#             bs_cap_flr_let = caplet_prices_black[i]
#             print("%d %9.6f %9.6f  %9.6f %9.6f %9.6f"% (i, sum_cap1_f[i] * 100.0,
#                                                          sum_flr1_f[i] * 100.0,
#                                                          sum_cap_nf[i] * 100.0,
#                                                          sum_flr_nf[i] * 100.0,
#                                                          bs_cap_flr_let * 100.0))

########################################################################################


########################################################################################


def test_HullBookExamples():
    """Examining examples on page 770 of Hull OFODS
    Last cap product has caplet starting in 10 years so we have to model
    the forward curve from time 0 out to 11 forwards, not 10 forwards.
    We have to model forward rates 0-1, 1-2, 2-3, ..., 10-11"""

    verbose = True

    # We go out 11 periods because last caplet resets in 10 years
    num_fwds = 11
    dt = 1.00
    taus = np.array([dt] * num_fwds)
    seed = 438

    r = 0.05127
    fwd0 = np.zeros(num_fwds)
    for i in range(0, num_fwds):
        fwd0[i] = r

    num_paths = 500000
    spread = 0.0025  # basis points

    test_cases.header("COMMENTS", "VALUES")

    ###########################################################################
    # HULL TABLE 32.1
    ###########################################################################

    use_sobol = 1
    numeraire_index = 0

    # We need the volatility for the forward rates out to the one starting in
    # 10 years. So we have 11 elements. The one starting today has zero vol.
    num_factors = 1
    gammas1FList = [
        0.00,
        0.1550,
        0.2063674,
        0.1720986,
        0.1721993,
        0.1524579,
        0.1414779,
        0.1297711,
        0.1381053,
        0.135955,
        0.1339842,
    ]
    gammas1F = np.array(gammas1FList)

    # One factor model
    fwds1F = lmm_simulate_fwds_1f(
        num_fwds,
        num_paths,
        numeraire_index,
        fwd0,
        gammas1F,
        taus,
        use_sobol,
        seed,
    )

    #    lmm_print_forwards(fwds1F)

    v_ratchet_caplets = (
        lmm_ratchet_caplet_pricer(spread, num_fwds, num_paths, fwd0, fwds1F, taus)
        * 100.0
    )

    hull_ratchet_caplets1_f = [
        0.00,
        0.196,
        0.207,
        0.201,
        0.194,
        0.187,
        0.1890,
        0.172,
        0.167,
        0.160,
        0.153,
    ]
    hull_ratchet_caplets1_f = np.array(hull_ratchet_caplets1_f)

    if verbose:
        test_cases.banner("Ratchet ONE FACTOR IMPLEMENTATION")
        test_cases.print("FINANCEPY GETS:", v_ratchet_caplets)
        test_cases.print("HULL GETS:", hull_ratchet_caplets1_f)

    check_vector_differences(v_ratchet_caplets, hull_ratchet_caplets1_f, 1e-2)

    v_sticky_caplets = (
        lmm_sticky_caplet_pricer(spread, num_fwds, num_paths, fwd0, fwds1F, taus) * 100.0
    )

    hull_sticky_caplets1_f = [
        0.0,
        0.196,
        0.336,
        0.412,
        0.458,
        0.484,
        0.498,
        0.502,
        0.501,
        0.497,
        0.488,
    ]

    if verbose:
        test_cases.banner("STICKY CAPLETS ONE FACTOR IMPLEMENTATION")
        test_cases.print("FINANCEPY GETS:", v_sticky_caplets)
        test_cases.print("HULL GETS:", hull_sticky_caplets1_f)

    check_vector_differences(v_sticky_caplets, hull_sticky_caplets1_f, 1e-2)

    num_factors = 1
    lambdas1FList = [
        [
            0.0,
            0.1550,
            0.2064,
            0.1721,
            0.1722,
            0.1525,
            0.1415,
            0.1298,
            0.1381,
            0.1360,
            0.1340,
        ]
    ]
    lambdas1F = np.array(lambdas1FList)

    # One factor model
    fwds_mf = lmm_simulate_fwds_mf(
        num_fwds,
        num_factors,
        num_paths,
        numeraire_index,
        fwd0,
        lambdas1F,
        taus,
        use_sobol,
        seed,
    )

    v_ratchet_caplets = (
        lmm_ratchet_caplet_pricer(spread, num_fwds, num_paths, fwd0, fwds_mf, taus)
        * 100.0
    )

    hull_ratchet_caplets1_f = [
        0.0,
        0.196,
        0.207,
        0.201,
        0.194,
        0.187,
        0.1890,
        0.172,
        0.167,
        0.160,
        0.153,
    ]

    if verbose:
        test_cases.banner("RATCHET - NUM FACTORS 1F")
        test_cases.print("FINANCEPY GETS:", v_ratchet_caplets)
        test_cases.print("HULL GETS:", hull_ratchet_caplets1_f)

    check_vector_differences(v_ratchet_caplets, hull_ratchet_caplets1_f, 1e-2)

    v_sticky_caplets = (
        lmm_sticky_caplet_pricer(spread, num_fwds, num_paths, fwd0, fwds_mf, taus) * 100.0
    )

    hull_sticky_caplets1_f = [
        0.00,
        0.196,
        0.336,
        0.412,
        0.458,
        0.484,
        0.498,
        0.502,
        0.501,
        0.497,
        0.488,
    ]

    check_vector_differences(v_sticky_caplets, hull_sticky_caplets1_f, 1e-2)

    if verbose:
        test_cases.banner("STICKY RATCHET - NUM FACTORS 1")
        test_cases.print("FINANCEPY GETS:", v_sticky_caplets)
        test_cases.print("HULL GETS:", hull_sticky_caplets1_f)

    num_factors = 2
    lambdas2FList = [
        [
            0.00,
            0.1410,
            0.1952,
            0.1678,
            0.1711,
            0.1525,
            0.1406,
            0.1265,
            0.1306,
            0.1236,
            0.1163,
        ],
        [
            0.00,
            -0.0645,
            -0.0670,
            -0.0384,
            -0.0196,
            0.00,
            0.0161,
            0.0289,
            0.0448,
            0.0565,
            0.0665,
        ],
    ]
    lambdas2F = np.array(lambdas2FList)

    # Two factor model
    fwds2F = lmm_simulate_fwds_mf(
        num_fwds,
        num_factors,
        num_paths,
        numeraire_index,
        fwd0,
        lambdas2F,
        taus,
        use_sobol,
        seed,
    )

    v_ratchet_caplets = (
        lmm_ratchet_caplet_pricer(spread, num_fwds, num_paths, fwd0, fwds2F, taus)
        * 100.0
    )
    hull_ratchet_caplets2_f = [
        0.00,
        0.194,
        0.207,
        0.205,
        0.198,
        0.193,
        0.189,
        0.180,
        0.174,
        0.168,
        0.162,
    ]

    if verbose:
        test_cases.banner("RATCHET - NUM FACTORS:2")
        test_cases.print("FINANCEPY GETS:", v_ratchet_caplets)
        test_cases.print("HULL GETS:", hull_ratchet_caplets2_f)

    check_vector_differences(v_ratchet_caplets, hull_ratchet_caplets2_f, 1e-2)

    v_sticky_caplets = (
        lmm_sticky_caplet_pricer(spread, num_fwds, num_paths, fwd0, fwds2F, taus) * 100.0
    )

    hull_sticky_caplets2_f = [
        0.00,
        0.196,
        0.334,
        0.413,
        0.462,
        0.492,
        0.512,
        0.520,
        0.523,
        0.523,
        0.519,
    ]

    if verbose:
        test_cases.banner("STICKY RATCHET - NUM FACTORS:2")
        test_cases.print("FINANCEPY GETS:", v_sticky_caplets)
        test_cases.print("HULL GETS:", hull_sticky_caplets2_f)

    check_vector_differences(v_sticky_caplets, hull_sticky_caplets2_f, 1e-2)

    num_factors = 3
    lambdas3FList = [
        [
            0.00,
            0.1365,
            0.1928,
            0.1672,
            0.1698,
            0.1485,
            0.1395,
            0.1261,
            0.1290,
            0.1197,
            0.1097,
        ],
        [
            0.0,
            -0.0662,
            -0.0702,
            -0.0406,
            -0.0206,
            0.00,
            0.0169,
            0.0306,
            0.0470,
            0.0581,
            0.0666,
        ],
        [
            0.0,
            0.0319,
            0.0225,
            0.000,
            -0.0198,
            -0.0347,
            -0.0163,
            0.000,
            0.0151,
            0.0280,
            0.0384,
        ],
    ]
    lambdas3F = np.array(lambdas3FList)

    # Three factor model
    fwds3F = lmm_simulate_fwds_mf(
        num_fwds,
        num_factors,
        num_paths,
        numeraire_index,
        fwd0,
        lambdas3F,
        taus,
        use_sobol,
        seed,
    )

    hull_ratchet_caplets3_f = [
        0.00,
        0.194,
        0.207,
        0.205,
        0.198,
        0.193,
        0.189,
        0.180,
        0.174,
        0.168,
        0.162,
    ]

    v_ratchet_caplets = (
        lmm_ratchet_caplet_pricer(spread, num_fwds, num_paths, fwd0, fwds3F, taus)
        * 100.0
    )

    if verbose:
        test_cases.banner("RATCHET - NUM FACTORS:3")
        test_cases.print("FINANCEPY GETS:", v_ratchet_caplets)
        test_cases.print("HULL GETS:", hull_ratchet_caplets3_f)

    check_vector_differences(v_ratchet_caplets, hull_ratchet_caplets3_f, 1e-2)

    v_sticky_caplets = (
        lmm_sticky_caplet_pricer(spread, num_fwds, num_paths, fwd0, fwds3F, taus) * 100.0
    )

    hull_sticky_caplets3_f = [
        0.00,
        0.195,
        0.336,
        0.418,
        0.472,
        0.506,
        0.524,
        0.533,
        0.537,
        0.537,
        0.534,
    ]

    if verbose:
        test_cases.banner("STICKY RATCHET - NUM FACTORS:3")
        test_cases.print("FINANCEPY GETS:", v_sticky_caplets)
        test_cases.print("HULL GETS:", hull_sticky_caplets3_f)

    check_vector_differences(v_sticky_caplets, hull_sticky_caplets3_f, 1e-2)


########################################################################################


""" def test_Swap():

    num_fwds = 40
    num_paths = 10000
    dt = 0.25
    taus = np.array([dt] * num_fwds)
    seed = 143

    r = 0.04
    fwd0 = get_forward_curve(num_fwds, r)

    fwd_rateVol = 0.40
    zetas = get_vol_curve(num_fwds, dt, fwd_rateVol)
    correl = get_correlation_matrix(num_fwds, 100.0)

    fwds = lmm_simulate_fwds_nf(num_fwds, num_paths, fwd0, zetas, correl, taus, seed)

    start = time.time()
    cpn = 0.05
    num_periods = 10
    swap = lmm_swap_pricer(cpn, num_periods, num_paths, fwd0, fwds, taus)
    end = time.time()
    print("PRICER Period:", end - start)

    print(swap) """

########################################################################################


########################################################################################


def fwdfwd_correlation(fwds):

    num_paths = len(fwds)
    num_fwds = len(fwds[0])
    start = time.time()
    fwd_corr = lmm_fwd_fwd_correlation(num_fwds, num_paths, 1, fwds)
    end = time.time()


#    print("CORR Period:", end - start)
#    print(fwd_corr)

########################################################################################


test_HullBookExamples()
# test_CapsFloors()
# test_Swaptions()
test_cases.compare_test_cases()

########################################################################################

