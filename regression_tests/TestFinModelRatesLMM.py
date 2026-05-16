# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

import time as time
import numpy as np

import add_fp_to_path

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

from FinTestCases import FinTestCases, global_test_case_mode

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


""" def getVolCurve(num_fwds, dt, flatVol=None):

    value_dt = Date(1, 1, 2020)

    capVolDates = []
    capletVolTenor = "1Y"
    num_periods = 10
    capletDt = value_dt

    capVolDates.append(value_dt)
    for _ in range(0, num_periods):
        capletDt = capletDt.add_tenor(capletVolTenor)
        capVolDates.append(capletDt)

    if flatVol is None:
        capVolatilities = [0.0, 15.50, 18.25, 17.91, 17.74, 17.27,
                           16.79, 16.30, 16.01, 15.76, 15.54]
        capVolatilities = np.array(capVolatilities)/100.0
    else:
        capVolatilities = [flatVol] * (num_periods+1)
        capVolatilities = np.array(capVolatilities)
        capVolatilities[0] = 0.0

    dc_type = DayCountTypes.ACT_ACT_ISDA
    vol_curve = IborCapVolCurve(value_dt,
                                   capVolDates,
                                   capVolatilities,
                                   dc_type)

    zetas = np.zeros(num_fwds)
    t = 0.0
    for ix in range(0, num_fwds):
        t = t + dt
        zetas[ix] = vol_curve.capFloorletVol(t)

    print(zetas)
    return zetas """

########################################################################################


def get_forward_curve(num_fwds, r):

    fwd0 = np.zeros(num_fwds)
    for i in range(0, num_fwds):
        fwd0[i] = r

    fwd0 = np.array(fwd0)
    return fwd0


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
#     zetas = getVolCurve(num_fwds, dt, fwd_rateVol)

#     seed = 1489
#     num_paths = 2000  # 100000
#     fwdsNF = LMMSimulateFwdsNF(num_fwds, num_paths, fwd0,
#                                zetas, correl, taus, seed)
#     strike = r
#     PAYSwaption = 1
#     use_sobol = 0
#     numeraire_index = 0

#     fwds1_f = LMMSimulateFwds1F(num_fwds, num_paths, numeraire_index, fwd0,
#                                zetas, taus, use_sobol, seed)

#     for iExp in range(1, 10):

#         t_exp = float(iExp)
#         a = int(2*t_exp)
#         print(a, b)

#         swaption_price1F = LMMSwaptionPricer(strike, a, b, num_paths,
#                                             fwd0, fwds1_f, taus, PAYSwaption)

#         swaption_priceNF = LMMSwaptionPricer(strike, a, b, num_paths,
#                                             fwd0, fwdsNF, taus, PAYSwaption)

#         swaption_vol = LMMSwaptionVolApprox(a, b, fwd0, taus, zetas, correl)

#         swapVolSim1F = LMMSimSwaptionVol(a, b, fwd0, fwds1_f, taus)
#         swapVolSimNF = LMMSimSwaptionVol(a, b, fwd0, fwdsNF, taus)

#         value_dt = Date(1, 1, 2010)
#         libor_curve = FinDiscountCurveFlat(value_dt, r,
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
#         swaptionType = IborSwaptionTypes.PAY
#         swaption = IborSwaption(settle_dt,
#                                     exercise_dt,
#                                     maturity_dt,
#                                     swaptionType,
#                                     fixed_cpn,
#                                     fixed_freq_type,
#                                     fixed_dc_type,
#                                     notional,
#                                     float_freq_type,
#                                     float_dc_type)

#         model = Black(swaption_vol)
#         blackSwaptionPrice = swaption.value(value_dt, libor_curve, model)

#         print("K:%6.5f t_exp:%8.2f FwdVol:%9.5f SimVol1F:%9.5f SimVolNF:%9.5f RebVol:%9.5f SimPx1F:%9.5f SimPxNF:%9.5f Black Px:%9.5f"
#               % (strike, t_exp, fwd_rateVol, swapVolSim1F, swapVolSimNF, swaption_vol,
#                  swaption_price1F, swaption_priceNF, blackSwaptionPrice))

# #        print(swaption)


# def test_CapsFloors():

#     num_fwds = 40
#     num_paths = 100000
#     dt = 0.25
#     taus = np.array([dt] * num_fwds)
#     seeds = [42] #, 143, 101, 993, 9988]

#     r = 0.04
#     fwd0 = get_forward_curve(num_fwds, r)

#     fwd_rateVol = 0.20
#     zetas = getVolCurve(num_fwds, dt, fwd_rateVol)

#     correl = get_correlation_matrix(num_fwds, 100.0, dt)

#     # At the money
#     K = r
#     capletPricesBlack = LMMPriceCapsBlack(fwd0, zetas, num_fwds, K, taus)

#     num_factors = 1
#     numeraire_index = 1
#     use_sobol = 1

#     # Examine variance for different seeds
#     for seed in seeds:

#         print("=============================================================")
#         print("Seed:", seed)

#         fwds1_f = LMMSimulateFwds1F(num_fwds, num_paths, numeraire_index, fwd0,
#                                    zetas, taus, use_sobol, seed)

#         sumCap1F = LMMCapFlrPricer(num_fwds, num_paths, K, fwd0, fwds1_f, taus, 1)
#         sumFlr1F = LMMCapFlrPricer(num_fwds, num_paths, K, fwd0, fwds1_f, taus, 0)

#         fwdsNF = LMMSimulateFwdsNF(num_fwds, num_paths, fwd0,
#                                    zetas, correl, taus, seed)

#         sumCapNF = LMMCapFlrPricer(num_fwds, num_paths, K, fwd0, fwdsNF, taus, 1)
#         sumFlrNF = LMMCapFlrPricer(num_fwds, num_paths, K, fwd0, fwdsNF, taus, 0)

#         print("i     CAPLET1F    FLRLET1F   CAPLETNF    FLRLET_NF   BS CAP")
#         for i in range(0, num_fwds):
#             bsCapFlrLet = capletPricesBlack[i]
#             print("%d %9.6f %9.6f  %9.6f %9.6f %9.6f"% (i, sumCap1F[i] * 100.0,
#                                                          sumFlr1F[i] * 100.0,
#                                                          sumCapNF[i] * 100.0,
#                                                          sumFlrNF[i] * 100.0,
#                                                          bsCapFlrLet * 100.0))

########################################################################################


def test_hull_book_examples():
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

    # HULL TABLE 32.1

    use_sobol = 1
    numeraire_index = 0

    # We need the volatility for the forward rates out to the one starting in
    # 10 years. So we have 11 elements. The one starting today has zero vol.
    num_factors = 1
    gammas1_f_list = [
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
    gammas1_f = np.array(gammas1_f_list)

    # One factor model
    fwds1_f = lmm_simulate_fwds_1f(
        num_fwds,
        num_paths,
        numeraire_index,
        fwd0,
        gammas1_f,
        taus,
        use_sobol,
        seed,
    )

    #    LMMPrintForwards(fwds1_f)

    v_ratchet_caplets = (
        lmm_ratchet_caplet_pricer(spread, num_fwds, num_paths, fwd0, fwds1_f, taus)
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
        lmm_sticky_caplet_pricer(spread, num_fwds, num_paths, fwd0, fwds1_f, taus)
        * 100.0
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
    lambdas1_f_list = [
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
    lambdas1_f = np.array(lambdas1_f_list)

    # One factor model
    fwds_mf = lmm_simulate_fwds_mf(
        num_fwds,
        num_factors,
        num_paths,
        numeraire_index,
        fwd0,
        lambdas1_f,
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
        lmm_sticky_caplet_pricer(spread, num_fwds, num_paths, fwd0, fwds_mf, taus)
        * 100.0
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
    lambdas2_f_list = [
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
    lambdas2_f = np.array(lambdas2_f_list)

    # Two factor model
    fwds2_f = lmm_simulate_fwds_mf(
        num_fwds,
        num_factors,
        num_paths,
        numeraire_index,
        fwd0,
        lambdas2_f,
        taus,
        use_sobol,
        seed,
    )

    v_ratchet_caplets = (
        lmm_ratchet_caplet_pricer(spread, num_fwds, num_paths, fwd0, fwds2_f, taus)
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
        lmm_sticky_caplet_pricer(spread, num_fwds, num_paths, fwd0, fwds2_f, taus)
        * 100.0
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
    lambdas3_f_list = [
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
    lambdas3_f = np.array(lambdas3_f_list)

    # Three factor model
    fwds3_f = lmm_simulate_fwds_mf(
        num_fwds,
        num_factors,
        num_paths,
        numeraire_index,
        fwd0,
        lambdas3_f,
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
        lmm_ratchet_caplet_pricer(spread, num_fwds, num_paths, fwd0, fwds3_f, taus)
        * 100.0
    )

    if verbose:
        test_cases.banner("RATCHET - NUM FACTORS:3")
        test_cases.print("FINANCEPY GETS:", v_ratchet_caplets)
        test_cases.print("HULL GETS:", hull_ratchet_caplets3_f)

    check_vector_differences(v_ratchet_caplets, hull_ratchet_caplets3_f, 1e-2)

    v_sticky_caplets = (
        lmm_sticky_caplet_pricer(spread, num_fwds, num_paths, fwd0, fwds3_f, taus)
        * 100.0
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


""" def test_Swap():

    num_fwds = 40
    num_paths = 10000
    dt = 0.25
    taus = np.array([dt] * num_fwds)
    seed = 143

    r = 0.04
    fwd0 = get_forward_curve(num_fwds, r)

    fwd_rateVol = 0.40
    zetas = getVolCurve(num_fwds, dt, fwd_rateVol)
    correl = get_correlation_matrix(num_fwds, 100.0)

    fwds = LMMSimulateFwdsNF(num_fwds, num_paths, fwd0, zetas, correl, taus, seed)

    start = time.time()
    cpn = 0.05
    num_periods = 10
    swap = LMMSwapPricer(cpn, num_periods, num_paths, fwd0, fwds, taus)
    end = time.time()
    print("PRICER Period:", end - start)

    print(swap) """

########################################################################################


def fwdfwd_correlation(fwds):

    num_paths = len(fwds)
    num_fwds = len(fwds[0])
    start = time.time()
    fwd_corr = lmm_fwd_fwd_correlation(num_fwds, num_paths, 1, fwds)
    end = time.time()


########################################################################################

#    print("CORR Period:", end - start)
#    print(fwd_corr)


test_hull_book_examples()
# test_CapsFloors()
# test_Swaptions()
test_cases.compare_test_cases()
