###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import numpy as np
from financepy.market.volatility.ibor_cap_vol_curve import IborCapVolCurve
from financepy.utils.date import Date
from financepy.utils.day_count import DayCountTypes
from financepy.models.black import Black
from financepy.market.curves.discount_curve_flat import DiscountCurveFlat
from financepy.utils.frequency import FrequencyTypes
from financepy.products.rates.ibor_swaption import SwapTypes
from financepy.products.rates.ibor_swaption import IborSwaption
from financepy.utils.global_types import FinCapFloorTypes
from financepy.products.rates.ibor_lmm_products import IborLMMProducts
from financepy.products.rates.ibor_cap_floor import IborCapFloor
from FinTestCases import FinTestCases, globalTestCaseMode
import sys
sys.path.append("..")


testCases = FinTestCases(__file__, globalTestCaseMode)

# This is in progress and needs to be completed

###############################################################################


# def test_Swaptions():
#     """ In progress and so not used. TODO. """

#     startYear = 2020
#     endYear = 2030
#     valuation_date = Date(1, 1, startYear)
#     exercise_date = Date(1, 1, 2023)
#     settlement_date = valuation_date
#     maturity_date = Date(1, 1, endYear)
#     fixed_coupon = 0.04

#     # DEFINE THE DISCOUNT CURVE
#     discount_curve = FinDiscountCurveFlat(valuation_date,
#                                          0.04,
#                                          FrequencyTypes.ANNUAL)

#     swaptionVol = 15.54

#     liborSwaption = IborSwaption(settlement_date,
#                                      exercise_date,
#                                      maturity_date,
#                                      IborSwaptionTypes.PAY,
#                                      fixed_coupon,
#                                      FrequencyTypes.ANNUAL,
#                                      DayCountTypes.ACT_360)

#     model = Black(swaptionVol/100.0)
#     v_BLK = liborSwaption.value(valuation_date, discount_curve, model)

#     dt = 0.5
#     texp = 3.0
#     tmat = 10.0
#     a = int(2*texp)
#     b = int(2*tmat)
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
#     useSobol = 0
#     numeraireIndex = 0

#     fwds1F = LMMSimulateFwds1F(numFwds, num_paths, numeraireIndex, fwd0,
#                                zetas, taus, useSobol, seed)

#     for iExp in range(1, 10):

#         texp = float(iExp)
#         a = int(2*texp)
#         print(a, b)

#         swaption_price1F = LMMSwaptionPricer(strike, a, b, num_paths,
#                                             fwd0, fwds1F, taus, PAYSwaption)

#         swaption_priceNF = LMMSwaptionPricer(strike, a, b, num_paths,
#                                             fwd0, fwdsNF, taus, PAYSwaption)

#         swaptionVol = LMMSwaptionVolApprox(a, b, fwd0, taus, zetas, correl)

#         swapVolSim1F = LMMSimSwaptionVol(a, b, fwd0, fwds1F, taus)
#         swapVolSimNF = LMMSimSwaptionVol(a, b, fwd0, fwdsNF, taus)

#         valuation_date = Date(1, 1, 2010)
#         libor_curve = FinDiscountCurveFlat(valuation_date, r,
#                                           FrequencyTypes.QUARTERLY)
#         settlement_date = valuation_date
#         exercise_date = settlement_date.add_months(a*3)
#         maturity_date = settlement_date.add_months(b*3)

#         fixed_coupon = strike
#         fixed_frequency_type = FrequencyTypes.QUARTERLY
#         fixed_day_count_type = DayCountTypes.ACT_ACT_ISDA
#         float_frequency_type = FrequencyTypes.QUARTERLY
#         float_day_count_type = DayCountTypes.ACT_ACT_ISDA
#         notional = 1.0

#         # Pricing a PAY
#         swaptionType = IborSwaptionTypes.PAY
#         swaption = IborSwaption(settlement_date,
#                                     exercise_date,
#                                     maturity_date,
#                                     swaptionType,
#                                     fixed_coupon,
#                                     fixed_frequency_type,
#                                     fixed_day_count_type,
#                                     notional,
#                                     float_frequency_type,
#                                     float_day_count_type)

#         model = Black(swaptionVol)
#         blackSwaptionPrice = swaption.value(valuation_date, libor_curve, model)

#         testCases.print("K:%6.5f texp:%8.2f FwdVol:%9.5f SimVol1F:%9.5f " +
#                         " SimVolNF:%9.5f RebVol:%9.5f SimPx1F:%9.5f SimPxNF:%9.5f Black Px:%9.5f"
#               % (strike, texp, fwd_rateVol, swapVolSim1F, swapVolSimNF,
#                  swaptionVol, swaption_price1F, swaption_priceNF,
#                  blackSwaptionPrice))

###############################################################################


# def test_CapsFloors():

#     # Define the CAP
#     # The maturity date is in 10 years so the last caplet start time is in 9
#     # years which in our convention means we are modelling 10 forwards
#     startYear = 2020
#     endYear = 2030
#     valuation_date = Date(1, 1, startYear)
#     settlement_date = valuation_date
#     capMaturityDate = Date(1, 1, endYear)
#     freq_type = FrequencyTypes.ANNUAL
#     day_count_type = DayCountTypes.ACT_360
#     capFloorRate = 0.04

#     # DEFINE THE DISCOUNT CURVE
#     discount_curve = FinDiscountCurveFlat(valuation_date,
#                                          0.04,
#                                          FrequencyTypes.ANNUAL)

#     capVol = 15.54

#     liborCap = IborCapFloor(settlement_date,
#                                 capMaturityDate,
#                                 IborCapFloorTypes.CAP,
#                                 capFloorRate,
#                                 None,
#                                 FrequencyTypes.ANNUAL,
#                                 DayCountTypes.ACT_360)

#     model = Black(capVol/100.0)
#     v_BLK = liborCap.value(valuation_date, discount_curve, model)

#     ###########################################################################
#     # LMM VALUATION
#     ###########################################################################

#     lmmProducts = IborLMMProducts(settlement_date,
#                                       capMaturityDate,
#                                       freq_type,
#                                       day_count_type)

#     # Set up forward rate vol structure
#     capVolDates = []
#     capletVolTenor = "1Y"
#     capletDt = valuation_date
#     numForwards = endYear - startYear

#     # Capvol dates has numForwards + 1 elements including today
#     capVolDates.append(valuation_date)
#     for i in range(0, numForwards):
#         capletDt = capletDt.add_tenor(capletVolTenor)
#         capVolDates.append(capletDt)

#     # Capvol dates has numForwards + 1 elements including zero today
#     capVolatilities = [capVol] * (numForwards+1)
#     capVolatilities[0] = 0.0
#     capVolatilities = np.array(capVolatilities)/100.0

#     day_count_type = DayCountTypes.ACT_ACT_ISDA
#     volCurve = IborCapVolCurve(valuation_date,
#                                    capVolDates,
#                                    capVolatilities,
#                                    day_count_type)

#     lambdas2FList = [[0.00, 0.1410, 0.1952, 0.1678, 0.1711, 0.1525,
#                       0.1406, 0.1265, 0.1306, 0.1236],
#                      [0.00, -0.0645, -0.0670, -0.0384, -0.0196, 0.00,
#                      0.0161, 0.0289, 0.0448, 0.0565]]
#     lambdas2F = np.array(lambdas2FList)

#     # Simulate paths of future Libor rates
#     numFactors = 1

#     testCases.header("NUMPATHS", "VLMM", "VBLK", "ERROR")

#     for num_paths in [10000, 20000, 50000, 100000, 200000, 400000, 1000000]:

#         if numFactors == 1:
#             lmmProducts.simulate1F(discount_curve, volCurve, num_paths, 0, True)
#         elif numFactors == 2:
#             lmmProducts.simulateMF(discount_curve, numFactors, lambdas2F,
#                                    num_paths, 0, True)

#         v_lmm = lmmProducts.valueCapFloor(settlement_date,
#                                           capMaturityDate,
#                                           IborCapFloorTypes.CAP,
#                                           capFloorRate,
#                                           FrequencyTypes.ANNUAL,
#                                           DayCountTypes.ACT_360)

#         err = v_lmm - v_BLK
#         testCases.print(num_paths, v_lmm, v_BLK, err)

###############################################################################


# test_CapsFloors()
# test_Swaptions()
testCases.compareTestCases()
