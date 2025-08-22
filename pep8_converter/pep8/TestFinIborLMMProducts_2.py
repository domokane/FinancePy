# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from FinTestCases import FinTestCases, global_test_case_mode
from financepy.products.rates.ibor_cap_floor import IborCapFloor
from financepy.products.rates.ibor_lmm_products import IborLMMProducts
from financepy.utils.global_types import FinCapFloorTypes
from financepy.products.rates.ibor_swaption import IborSwaption
from financepy.products.rates.ibor_swaption import SwapTypes
from financepy.utils.frequency import FrequencyTypes
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.models.black import Black
from financepy.utils.day_count import DayCountTypes
from financepy.utils.date import Date
from financepy.market.volatility.ibor_cap_vol_curve import IborCapVolCurve
import numpy as np

import sys

sys.path.append("..")


test_cases = FinTestCases(__file__, global_test_case_mode)

# This is in progress and needs to be completed



# def test_Swaptions():
#     """ In progress and so not used. TODO. """

#     start_year = 2020
#     end_year = 2030
#     value_dt = Date(1, 1, start_year)
#     exercise_dt = Date(1, 1, 2023)
#     settle_dt = value_dt
#     maturity_dt = Date(1, 1, end_year)
#     fixed_cpn = 0.04

#     # DEFINE THE DISCOUNT CURVE
#     discount_curve = fin_discount_curve_flat(value_dt,
#                                          0.04,
#                                          FrequencyTypes.ANNUAL)

#     swaption_vol = 15.54

#     libor_swaption = IborSwaption(settle_dt,
#                                      exercise_dt,
#                                      maturity_dt,
#                                      ibor_swaption_types.PAY,
#                                      fixed_cpn,
#                                      FrequencyTypes.ANNUAL,
#                                      DayCountTypes.ACT_360)

#     model = Black(swaption_vol/100.0)
#     v_BLK = libor_swaption.value(value_dt, discount_curve, model)

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
#     num_paths = 2000 # 100000
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

#         test_cases.print("K:%6.5f t_exp:%8.2f fwd_vol:%9.5f sim_vol1_f:%9.5f " +
#                         " sim_vol_nf:%9.5f reb_vol:%9.5f sim_px1_f:%9.5f sim_px_nf:%9.5f Black Px:%9.5f"
#               % (strike, t_exp, fwd_rateVol, swap_vol_sim1_f, swap_vol_sim_nf,
#                  swaption_vol, swaption_price1F, swaption_priceNF,
#                  black_swaption_price))



# def test_CapsFloors():

#     # Define the CAP
#     # The maturity date is in 10 years so the last caplet start time is in 9
#     # years which in our convention means we are modelling 10 forwards
#     start_year = 2020
#     end_year = 2030
#     value_dt = Date(1, 1, start_year)
#     settle_dt = value_dt
#     cap_maturity_date = Date(1, 1, end_year)
#     freq_type = FrequencyTypes.ANNUAL
#     dc_type = DayCountTypes.ACT_360
#     cap_floor_rate = 0.04

#     # DEFINE THE DISCOUNT CURVE
#     discount_curve = fin_discount_curve_flat(value_dt,
#                                          0.04,
#                                          FrequencyTypes.ANNUAL)

#     cap_vol = 15.54

#     libor_cap = IborCapFloor(settle_dt,
#                                 cap_maturity_date,
#                                 ibor_cap_floor_types.CAP,
#                                 cap_floor_rate,
#                                 None,
#                                 FrequencyTypes.ANNUAL,
#                                 DayCountTypes.ACT_360)

#     model = Black(cap_vol/100.0)
#     v_BLK = libor_cap.value(value_dt, discount_curve, model)

#     ###########################################################################
#     # LMM VALUATION
#     ###########################################################################

#     lmm_products = IborLMMProducts(settle_dt,
#                                       cap_maturity_date,
#                                       freq_type,
#                                       dc_type)

#     # Set up forward rate vol structure
#     cap_vol_dates = []
#     caplet_vol_tenor = "1Y"
#     caplet_dt = value_dt
#     num_fwds = end_year - start_year

#     # Capvol dates has num_fwds + 1 elements including today
#     cap_vol_dates.append(value_dt)
#     for i in range(0, num_fwds):
#         caplet_dt = caplet_dt.add_tenor(caplet_vol_tenor)
#         cap_vol_dates.append(caplet_dt)

#     # Capvol dates has num_fwds + 1 elements including zero today
#     cap_volatilities = [cap_vol] * (num_fwds+1)
#     cap_volatilities[0] = 0.0
#     cap_volatilities = np.array(cap_volatilities)/100.0

#     dc_type = DayCountTypes.ACT_ACT_ISDA
#     vol_curve = IborCapVolCurve(value_dt,
#                                    cap_vol_dates,
#                                    cap_volatilities,
#                                    dc_type)

#     lambdas2FList = [[0.00, 0.1410, 0.1952, 0.1678, 0.1711, 0.1525,
#                       0.1406, 0.1265, 0.1306, 0.1236],
#                      [0.00, -0.0645, -0.0670, -0.0384, -0.0196, 0.00,
#                      0.0161, 0.0289, 0.0448, 0.0565]]
#     lambdas2F = np.array(lambdas2FList)

#     # Simulate paths of future Libor rates
#     num_factors = 1

#     test_cases.header("NUMPATHS", "VLMM", "VBLK", "ERROR")

#     for num_paths in [10000, 20000, 50000, 100000, 200000, 400000, 1000000]:

#         if num_factors == 1:
#             lmm_products.simulate1F(discount_curve, vol_curve, num_paths, 0, True)
#         elif num_factors == 2:
#             lmm_products.simulate_mf(discount_curve, num_factors, lambdas2F,
#                                    num_paths, 0, True)

#         v_lmm = lmm_products.value_cap_floor(settle_dt,
#                                           cap_maturity_date,
#                                           ibor_cap_floor_types.CAP,
#                                           cap_floor_rate,
#                                           FrequencyTypes.ANNUAL,
#                                           DayCountTypes.ACT_360)

#         err = v_lmm - v_BLK
#         test_cases.print(num_paths, v_lmm, v_BLK, err)



# test_CapsFloors()
# test_Swaptions()
test_cases.compare_test_cases()
