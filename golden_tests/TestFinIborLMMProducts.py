###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from FinTestCases import FinTestCases, globalTestCaseMode
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


test_cases = FinTestCases(__file__, globalTestCaseMode)

# This is in progress and needs to be completed

###############################################################################


# def test_Swaptions():
#     """ In progress and so not used. TODO. """

#     startYear = 2020
#     endYear = 2030
#     value_dt = Date(1, 1, startYear)
#     exercise_dt = Date(1, 1, 2023)
#     settle_dt = value_dt
#     maturity_dt = Date(1, 1, endYear)
#     fixed_cpn = 0.04

#     # DEFINE THE DISCOUNT CURVE
#     discount_curve = FinDiscountCurveFlat(value_dt,
#                                          0.04,
#                                          FrequencyTypes.ANNUAL)

#     swaption_vol = 15.54

#     liborSwaption = IborSwaption(settle_dt,
#                                      exercise_dt,
#                                      maturity_dt,
#                                      IborSwaptionTypes.PAY,
#                                      fixed_cpn,
#                                      FrequencyTypes.ANNUAL,
#                                      DayCountTypes.ACT_360)

#     model = Black(swaption_vol/100.0)
#     v_BLK = liborSwaption.value(value_dt, discount_curve, model)

#     dt = 0.5
#     t_exp = 3.0
#     t_mat = 10.0
#     a = int(2*t_exp)
#     b = int(2*t_mat)
#     numFwds = 20
#     taus = np.array([dt] * numFwds)

#     r = 0.05
#     fwd0 = getForwardCurve(numFwds, r)
#     correl = getCorrelationMatrix(numFwds, 0.00000000000001, dt)

#     fwd_rateVol = 0.20
#     zetas = getVolCurve(numFwds, dt, fwd_rateVol)

#     seed = 1489
#     num_paths = 2000 # 100000
#     fwdsNF = LMMSimulateFwdsNF(numFwds, num_paths, fwd0,
#                                zetas, correl, taus, seed)
#     strike = r
#     PAYSwaption = 1
#     use_sobol = 0
#     numeraire_index = 0

#     fwds1F = LMMSimulateFwds1F(numFwds, num_paths, numeraire_index, fwd0,
#                                zetas, taus, use_sobol, seed)

#     for iExp in range(1, 10):

#         t_exp = float(iExp)
#         a = int(2*t_exp)
#         print(a, b)

#         swaption_price1F = LMMSwaptionPricer(strike, a, b, num_paths,
#                                             fwd0, fwds1F, taus, PAYSwaption)

#         swaption_priceNF = LMMSwaptionPricer(strike, a, b, num_paths,
#                                             fwd0, fwdsNF, taus, PAYSwaption)

#         swaption_vol = LMMSwaptionVolApprox(a, b, fwd0, taus, zetas, correl)

#         swapVolSim1F = LMMSimSwaptionVol(a, b, fwd0, fwds1F, taus)
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

#         test_cases.print("K:%6.5f t_exp:%8.2f FwdVol:%9.5f SimVol1F:%9.5f " +
#                         " SimVolNF:%9.5f RebVol:%9.5f SimPx1F:%9.5f SimPxNF:%9.5f Black Px:%9.5f"
#               % (strike, t_exp, fwd_rateVol, swapVolSim1F, swapVolSimNF,
#                  swaption_vol, swaption_price1F, swaption_priceNF,
#                  blackSwaptionPrice))

###############################################################################


# def test_CapsFloors():

#     # Define the CAP
#     # The maturity date is in 10 years so the last caplet start time is in 9
#     # years which in our convention means we are modelling 10 forwards
#     startYear = 2020
#     endYear = 2030
#     value_dt = Date(1, 1, startYear)
#     settle_dt = value_dt
#     capMaturityDate = Date(1, 1, endYear)
#     freq_type = FrequencyTypes.ANNUAL
#     dc_type = DayCountTypes.ACT_360
#     capFloorRate = 0.04

#     # DEFINE THE DISCOUNT CURVE
#     discount_curve = FinDiscountCurveFlat(value_dt,
#                                          0.04,
#                                          FrequencyTypes.ANNUAL)

#     capVol = 15.54

#     liborCap = IborCapFloor(settle_dt,
#                                 capMaturityDate,
#                                 IborCapFloorTypes.CAP,
#                                 capFloorRate,
#                                 None,
#                                 FrequencyTypes.ANNUAL,
#                                 DayCountTypes.ACT_360)

#     model = Black(capVol/100.0)
#     v_BLK = liborCap.value(value_dt, discount_curve, model)

#     ###########################################################################
#     # LMM VALUATION
#     ###########################################################################

#     lmmProducts = IborLMMProducts(settle_dt,
#                                       capMaturityDate,
#                                       freq_type,
#                                       dc_type)

#     # Set up forward rate vol structure
#     capVolDates = []
#     capletVolTenor = "1Y"
#     capletDt = value_dt
#     num_fwds = endYear - startYear

#     # Capvol dates has num_fwds + 1 elements including today
#     capVolDates.append(value_dt)
#     for i in range(0, num_fwds):
#         capletDt = capletDt.add_tenor(capletVolTenor)
#         capVolDates.append(capletDt)

#     # Capvol dates has num_fwds + 1 elements including zero today
#     capVolatilities = [capVol] * (num_fwds+1)
#     capVolatilities[0] = 0.0
#     capVolatilities = np.array(capVolatilities)/100.0

#     dc_type = DayCountTypes.ACT_ACT_ISDA
#     vol_curve = IborCapVolCurve(value_dt,
#                                    capVolDates,
#                                    capVolatilities,
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
#             lmmProducts.simulate1F(discount_curve, vol_curve, num_paths, 0, True)
#         elif num_factors == 2:
#             lmmProducts.simulateMF(discount_curve, num_factors, lambdas2F,
#                                    num_paths, 0, True)

#         v_lmm = lmmProducts.valueCapFloor(settle_dt,
#                                           capMaturityDate,
#                                           IborCapFloorTypes.CAP,
#                                           capFloorRate,
#                                           FrequencyTypes.ANNUAL,
#                                           DayCountTypes.ACT_360)

#         err = v_lmm - v_BLK
#         test_cases.print(num_paths, v_lmm, v_BLK, err)

###############################################################################


# test_CapsFloors()
# test_Swaptions()
test_cases.compareTestCases()
