##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import time as time
import numpy as np

import sys
sys.path.append("..")

from financepy.market.volatility.FinIborCapVolCurve import FinIborCapVolCurve
from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.models.FinModelBlack import FinModelBlack
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.products.rates.FinIborSwaption import FinSwapTypes
from financepy.products.rates.FinIborSwaption import FinIborSwaption
from financepy.finutils.FinHelperFunctions import checkVectorDifferences
from financepy.models.FinModelRatesLMM import LMMSimulateFwdsNF
from financepy.models.FinModelRatesLMM import LMMSimulateFwds1F
from financepy.models.FinModelRatesLMM import LMMSimulateFwdsMF
from financepy.models.FinModelRatesLMM import LMMSwaptionPricer
from financepy.models.FinModelRatesLMM import LMMSimSwaptionVol
from financepy.models.FinModelRatesLMM import LMMSwaptionVolApprox
from financepy.models.FinModelRatesLMM import LMMCapFlrPricer
from financepy.models.FinModelRatesLMM import LMMPriceCapsBlack
from financepy.models.FinModelRatesLMM import LMMSwapPricer
from financepy.models.FinModelRatesLMM import LMMFwdFwdCorrelation
from financepy.models.FinModelRatesLMM import LMMRatchetCapletPricer
from financepy.models.FinModelRatesLMM import LMMStickyCapletPricer

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def getCorrelationMatrix(numFwds, beta, dt):
    correl = np.zeros((numFwds, numFwds))  # indexes from 1 to p-1
    for i in range(0, numFwds):
        correl[i][i] = 1.0
        for j in range(i+1, numFwds):
            correl[i][j] = np.exp(-beta * (j - i) * dt)
            correl[j][i] = correl[i][j]

    return correl

###############################################################################


""" def getVolCurve(numFwds, dt, flatVol=None):

    valuationDate = FinDate(1, 1, 2020)

    capVolDates = []
    capletVolTenor = "1Y"
    numPeriods = 10
    capletDt = valuationDate

    capVolDates.append(valuationDate)
    for _ in range(0, numPeriods):
        capletDt = capletDt.addTenor(capletVolTenor)
        capVolDates.append(capletDt)

    if flatVol is None:
        capVolatilities = [0.0, 15.50, 18.25, 17.91, 17.74, 17.27,
                           16.79, 16.30, 16.01, 15.76, 15.54]
        capVolatilities = np.array(capVolatilities)/100.0
    else:
        capVolatilities = [flatVol] * (numPeriods+1)
        capVolatilities = np.array(capVolatilities)
        capVolatilities[0] = 0.0

    dayCountType = FinDayCountTypes.ACT_ACT_ISDA
    volCurve = FinIborCapVolCurve(valuationDate,
                                   capVolDates,
                                   capVolatilities,
                                   dayCountType)

    zetas = np.zeros(numFwds)
    t = 0.0
    for ix in range(0, numFwds):
        t = t + dt
        zetas[ix] = volCurve.capFloorletVol(t)

    print(zetas)
    return zetas """

###############################################################################


def getForwardCurve(numFwds, r):

    fwd0 = np.zeros(numFwds)
    for i in range(0, numFwds):
        fwd0[i] = r

    fwd0 = np.array(fwd0)
    return fwd0

###############################################################################


# def test_Swaptions():

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

#     fwdRateVol = 0.20
#     zetas = getVolCurve(numFwds, dt, fwdRateVol)

#     seed = 1489
#     numPaths = 2000  # 100000
#     fwdsNF = LMMSimulateFwdsNF(numFwds, numPaths, fwd0,
#                                zetas, correl, taus, seed)
#     strike = r
#     PAYSwaption = 1
#     useSobol = 0
#     numeraireIndex = 0

#     fwds1F = LMMSimulateFwds1F(numFwds, numPaths, numeraireIndex, fwd0,
#                                zetas, taus, useSobol, seed)

#     for iExp in range(1, 10):

#         texp = float(iExp)
#         a = int(2*texp)
#         print(a, b)

#         swaptionPrice1F = LMMSwaptionPricer(strike, a, b, numPaths,
#                                             fwd0, fwds1F, taus, PAYSwaption)

#         swaptionPriceNF = LMMSwaptionPricer(strike, a, b, numPaths,
#                                             fwd0, fwdsNF, taus, PAYSwaption)

#         swaptionVol = LMMSwaptionVolApprox(a, b, fwd0, taus, zetas, correl)

#         swapVolSim1F = LMMSimSwaptionVol(a, b, fwd0, fwds1F, taus)
#         swapVolSimNF = LMMSimSwaptionVol(a, b, fwd0, fwdsNF, taus)

#         valuationDate = FinDate(1, 1, 2010)
#         liborCurve = FinDiscountCurveFlat(valuationDate, r,
#                                           FinFrequencyTypes.QUARTERLY)

#         settlementDate = valuationDate
#         exerciseDate = settlementDate.addMonths(a*3)
#         maturityDate = settlementDate.addMonths(b*3)

#         fixedCoupon = strike
#         fixedFrequencyType = FinFrequencyTypes.QUARTERLY
#         fixedDayCountType = FinDayCountTypes.ACT_ACT_ISDA
#         floatFrequencyType = FinFrequencyTypes.QUARTERLY
#         floatDayCountType = FinDayCountTypes.ACT_ACT_ISDA
#         notional = 1.0

#         # Pricing a PAY
#         swaptionType = FinIborSwaptionTypes.PAY
#         swaption = FinIborSwaption(settlementDate,
#                                     exerciseDate,
#                                     maturityDate,
#                                     swaptionType,
#                                     fixedCoupon,
#                                     fixedFrequencyType,
#                                     fixedDayCountType,
#                                     notional,
#                                     floatFrequencyType,
#                                     floatDayCountType)

#         model = FinModelBlack(swaptionVol)
#         blackSwaptionPrice = swaption.value(valuationDate, liborCurve, model)

#         print("K:%6.5f texp:%8.2f FwdVol:%9.5f SimVol1F:%9.5f SimVolNF:%9.5f RebVol:%9.5f SimPx1F:%9.5f SimPxNF:%9.5f Black Px:%9.5f" 
#               % (strike, texp, fwdRateVol, swapVolSim1F, swapVolSimNF, swaptionVol,
#                  swaptionPrice1F, swaptionPriceNF, blackSwaptionPrice))

# #        print(swaption)

###############################################################################


# def test_CapsFloors():

#     numFwds = 40
#     numPaths = 100000
#     dt = 0.25
#     taus = np.array([dt] * numFwds)
#     seeds = [42] #, 143, 101, 993, 9988]

#     r = 0.04
#     fwd0 = getForwardCurve(numFwds, r)

#     fwdRateVol = 0.20
#     zetas = getVolCurve(numFwds, dt, fwdRateVol)

#     correl = getCorrelationMatrix(numFwds, 100.0, dt)

#     # At the money
#     K = r
#     capletPricesBlack = LMMPriceCapsBlack(fwd0, zetas, numFwds, K, taus)

#     numFactors = 1
#     numeraireIndex = 1
#     useSobol = 1

#     # Examine variance for different seeds
#     for seed in seeds:

#         print("=============================================================")
#         print("Seed:", seed)

#         fwds1F = LMMSimulateFwds1F(numFwds, numPaths, numeraireIndex, fwd0,
#                                    zetas, taus, useSobol, seed)

#         sumCap1F = LMMCapFlrPricer(numFwds, numPaths, K, fwd0, fwds1F, taus, 1)
#         sumFlr1F = LMMCapFlrPricer(numFwds, numPaths, K, fwd0, fwds1F, taus, 0)

#         fwdsNF = LMMSimulateFwdsNF(numFwds, numPaths, fwd0,
#                                    zetas, correl, taus, seed)

#         sumCapNF = LMMCapFlrPricer(numFwds, numPaths, K, fwd0, fwdsNF, taus, 1)
#         sumFlrNF = LMMCapFlrPricer(numFwds, numPaths, K, fwd0, fwdsNF, taus, 0)

#         print("i     CAPLET1F    FLRLET1F   CAPLETNF    FLRLET_NF   BS CAP")
#         for i in range(0, numFwds):
#             bsCapFlrLet = capletPricesBlack[i]
#             print("%d %9.6f %9.6f  %9.6f %9.6f %9.6f"% (i, sumCap1F[i] * 100.0,
#                                                          sumFlr1F[i] * 100.0,
#                                                          sumCapNF[i] * 100.0,
#                                                          sumFlrNF[i] * 100.0,
#                                                          bsCapFlrLet * 100.0))

###############################################################################


def test_HullBookExamples():
    ''' Examining examples on page 770 of Hull OFODS
        Last cap product has caplet starting in 10 years so we have to model
        the forward curve from time 0 out to 11 forwards, not 10 forwards.
        We have to model forward rates 0-1, 1-2, 2-3, ..., 10-11 '''

    verbose = True

    # We go out 11 periods because last caplet resets in 10 years
    numFwds = 11
    dt = 1.00
    taus = np.array([dt] * numFwds)
    seed = 438

    r = 0.05127
    fwd0 = np.zeros(numFwds)
    for i in range(0, numFwds):
        fwd0[i] = r

    numPaths = 500000
    spread = 0.0025  # basis points

    testCases.header("COMMENTS", "VALUES")

    ###########################################################################
    # HULL TABLE 32.1
    ###########################################################################

    useSobol = 1
    numeraireIndex = 0

    # We need the volatility for the forward rates out to the one starting in
    # 10 years. So we have 11 elements. The one starting today has zero vol.
    numFactors = 1
    gammas1FList = [0.00, 0.1550, 0.2063674, 0.1720986, 0.1721993, 0.1524579,
                    0.1414779, 0.1297711, 0.1381053, 0.135955, 0.1339842]
    gammas1F = np.array(gammas1FList)

    # One factor model
    fwds1F = LMMSimulateFwds1F(numFwds, numPaths, numeraireIndex, fwd0,
                               gammas1F, taus, useSobol, seed)

#    LMMPrintForwards(fwds1F)

    vRatchetCaplets = LMMRatchetCapletPricer(spread, numFwds, numPaths,
                                             fwd0, fwds1F, taus) * 100.0

    hullRatchetCaplets1F = [0.00, 0.196, 0.207, 0.201, 0.194, 0.187,
                            0.1890, 0.172, 0.167, 0.160, 0.153]
    hullRatchetCaplets1F = np.array(hullRatchetCaplets1F)

    if verbose:
        testCases.banner("Ratchet ONE FACTOR IMPLEMENTATION")
        testCases.print("FINANCEPY GETS:", vRatchetCaplets)
        testCases.print("HULL GETS:", hullRatchetCaplets1F)

    checkVectorDifferences(vRatchetCaplets, hullRatchetCaplets1F, 1e-2)

    vStickyCaplets = LMMStickyCapletPricer(spread, numFwds, numPaths,
                                           fwd0, fwds1F, taus) * 100.0

    hullStickyCaplets1F = [0.0, 0.196, 0.336, 0.412, 0.458, 0.484,
                           0.498, 0.502, 0.501, 0.497, 0.488]

    if verbose:
        testCases.banner("STICKY CAPLETS ONE FACTOR IMPLEMENTATION")
        testCases.print("FINANCEPY GETS:", vStickyCaplets)
        testCases.print("HULL GETS:", hullStickyCaplets1F)

    checkVectorDifferences(vStickyCaplets, hullStickyCaplets1F, 1e-2)

    numFactors = 1
    lambdas1FList = [[0.0, 0.1550, 0.2064, 0.1721, 0.1722, 0.1525,
                      0.1415, 0.1298, 0.1381, 0.1360, 0.1340]]
    lambdas1F = np.array(lambdas1FList)

    # One factor model
    fwdsMF = LMMSimulateFwdsMF(numFwds, numFactors, numPaths, numeraireIndex,
                               fwd0, lambdas1F, taus, useSobol, seed)

    vRatchetCaplets = LMMRatchetCapletPricer(spread, numFwds, numPaths,
                                             fwd0, fwdsMF, taus) * 100.0

    hullRatchetCaplets1F = [0.0, 0.196, 0.207, 0.201, 0.194, 0.187,
                            0.1890, 0.172, 0.167, 0.160, 0.153]

    if verbose:
        testCases.banner("RATCHET - NUM FACTORS 1F")
        testCases.print("FINANCEPY GETS:", vRatchetCaplets)
        testCases.print("HULL GETS:", hullRatchetCaplets1F)

    checkVectorDifferences(vRatchetCaplets, hullRatchetCaplets1F, 1e-2)

    vStickyCaplets = LMMStickyCapletPricer(spread, numFwds, numPaths,
                                           fwd0, fwdsMF, taus) * 100.0

    hullStickyCaplets1F = [0.00, 0.196, 0.336, 0.412, 0.458, 0.484, 0.498,
                           0.502, 0.501, 0.497, 0.488]

    checkVectorDifferences(vStickyCaplets, hullStickyCaplets1F, 1e-2)

    if verbose:
        testCases.banner("STICKY RATCHET - NUM FACTORS 1")
        testCases.print("FINANCEPY GETS:", vStickyCaplets)
        testCases.print("HULL GETS:", hullStickyCaplets1F)

    numFactors = 2
    lambdas2FList = [[0.00, 0.1410, 0.1952, 0.1678, 0.1711, 0.1525,
                      0.1406, 0.1265, 0.1306, 0.1236, 0.1163],
                     [0.00, -0.0645, -0.0670, -0.0384, -0.0196, 0.00,
                     0.0161, 0.0289, 0.0448, 0.0565, 0.0665]]
    lambdas2F = np.array(lambdas2FList)

    # Two factor model
    fwds2F = LMMSimulateFwdsMF(numFwds, numFactors, numPaths, numeraireIndex,
                               fwd0, lambdas2F, taus, useSobol, seed)

    vRatchetCaplets = LMMRatchetCapletPricer(spread, numFwds, numPaths,
                                             fwd0, fwds2F, taus) * 100.0
    hullRatchetCaplets2F = [0.00, 0.194, 0.207, 0.205, 0.198, 0.193,
                            0.189, 0.180, 0.174, 0.168, 0.162]

    if verbose:
        testCases.banner("RATCHET - NUM FACTORS:2")
        testCases.print("FINANCEPY GETS:", vRatchetCaplets)
        testCases.print("HULL GETS:", hullRatchetCaplets2F)

    checkVectorDifferences(vRatchetCaplets, hullRatchetCaplets2F, 1e-2)

    vStickyCaplets = LMMStickyCapletPricer(spread, numFwds, numPaths,
                                           fwd0, fwds2F, taus) * 100.0

    hullStickyCaplets2F = [0.00, 0.196, 0.334, 0.413, 0.462, 0.492,
                           0.512, 0.520, 0.523, 0.523, 0.519]

    if verbose:
        testCases.banner("STICKY RATCHET - NUM FACTORS:2")
        testCases.print("FINANCEPY GETS:", vStickyCaplets)
        testCases.print("HULL GETS:", hullStickyCaplets2F)

    checkVectorDifferences(vStickyCaplets, hullStickyCaplets2F, 1e-2)

    numFactors = 3
    lambdas3FList = [[0.00, 0.1365, 0.1928, 0.1672, 0.1698, 0.1485,
                     0.1395, 0.1261, 0.1290, 0.1197, 0.1097],
                     [0.0, -0.0662, -0.0702, -0.0406, -0.0206, 0.00,
                     0.0169, 0.0306, 0.0470, 0.0581, 0.0666],
                     [0.0, 0.0319, 0.0225, 0.000, -0.0198, -0.0347,
                     -0.0163, 0.000, 0.0151, 0.0280, 0.0384]]
    lambdas3F = np.array(lambdas3FList)

    # Three factor model
    fwds3F = LMMSimulateFwdsMF(numFwds, numFactors, numPaths, numeraireIndex,
                               fwd0, lambdas3F, taus, useSobol, seed)

    hullRatchetCaplets3F = [0.00, 0.194, 0.207, 0.205, 0.198, 0.193,
                            0.189, 0.180, 0.174, 0.168, 0.162]

    vRatchetCaplets = LMMRatchetCapletPricer(spread, numFwds, numPaths,
                                             fwd0, fwds3F, taus) * 100.0

    if verbose:
        testCases.banner("RATCHET - NUM FACTORS:3")
        testCases.print("FINANCEPY GETS:", vRatchetCaplets)
        testCases.print("HULL GETS:", hullRatchetCaplets3F)

    checkVectorDifferences(vRatchetCaplets, hullRatchetCaplets3F, 1e-2)

    vStickyCaplets = LMMStickyCapletPricer(spread, numFwds, numPaths,
                                           fwd0, fwds3F, taus) * 100.0

    hullStickyCaplets3F = [0.00, 0.195, 0.336, 0.418, 0.472, 0.506,
                           0.524, 0.533, 0.537, 0.537, 0.534]

    if verbose:
        testCases.banner("STICKY RATCHET - NUM FACTORS:3")
        testCases.print("FINANCEPY GETS:", vStickyCaplets)
        testCases.print("HULL GETS:", hullStickyCaplets3F)

    checkVectorDifferences(vStickyCaplets, hullStickyCaplets3F, 1e-2)

###############################################################################


""" def test_Swap():

    numFwds = 40
    numPaths = 10000
    dt = 0.25
    taus = np.array([dt] * numFwds)
    seed = 143

    r = 0.04
    fwd0 = getForwardCurve(numFwds, r)

    fwdRateVol = 0.40
    zetas = getVolCurve(numFwds, dt, fwdRateVol)
    correl = getCorrelationMatrix(numFwds, 100.0)

    fwds = LMMSimulateFwdsNF(numFwds, numPaths, fwd0, zetas, correl, taus, seed)

    start = time.time()
    cpn = 0.05
    numPeriods = 10
    swap = LMMSwapPricer(cpn, numPeriods, numPaths, fwd0, fwds, taus)
    end = time.time()
    print("PRICER Period:", end - start)

    print(swap) """

###############################################################################


def fwdfwdCorrelation(fwds):

    numPaths = len(fwds)
    numFwds = len(fwds[0])
    start = time.time()
    fwdCorr = LMMFwdFwdCorrelation(numFwds, numPaths, 1, fwds)
    end = time.time()
#    print("CORR Period:", end - start)
#    print(fwdCorr)

###############################################################################

test_HullBookExamples()
# test_CapsFloors()
# test_Swaptions()
testCases.compareTestCases()
