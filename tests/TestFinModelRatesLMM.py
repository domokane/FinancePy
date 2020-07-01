# -*- coding: utf-8 -*-

import time as time
from financepy.market.volatility.FinCapVolCurve import FinCapVolCurve
from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.models.FinModelBlack import FinModelBlack
from financepy.market.curves.FinFlatCurve import FinFlatCurve
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.products.libor.FinLiborSwaption import FinLiborSwaptionTypes
from financepy.products.libor.FinLiborSwaption import FinLiborSwaption
from financepy.finutils.FinHelperFunctions import checkVectorDifferences

from financepy.models.FinModelRatesLMM import LMMSimulateFwdsNF
from financepy.models.FinModelRatesLMM import LMMSimulateFwds1F
from financepy.models.FinModelRatesLMM import LMMSimulateFwdsMF
from financepy.models.FinModelRatesLMM import LMMSwaptionPricer
from financepy.models.FinModelRatesLMM import LMMSimSwaptionVol
from financepy.models.FinModelRatesLMM import LMMSwaptionVolApprox
from financepy.models.FinModelRatesLMM import LMMCapFlrPricer
from financepy.models.FinModelRatesLMM import priceCapsBlack
from financepy.models.FinModelRatesLMM import LMMSwapPricer
from financepy.models.FinModelRatesLMM import LMMFwdFwdCorrelation
from financepy.models.FinModelRatesLMM import LMMRatchetCapletPricer
from financepy.models.FinModelRatesLMM import LMMStickyCapletPricer
from financepy.models.FinModelRatesLMM import LMMFlexiCapPricer
from financepy.models.FinModelRatesLMM import LMMPrintForwards

import numpy as np

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


def getVolCurve(numFwds, dt, flatVol=None):

    valuationDate = FinDate(1, 1, 2020)

    capVolDates = []
    capletVolTenor = "1Y"
    numPeriods = 10
    capletDt = valuationDate

    for i in range(0, numPeriods):
        capletDt = capletDt.addTenor(capletVolTenor)
        capVolDates.append(capletDt)

    if flatVol is None:
        capVolatilities = [15.50, 18.25, 17.91, 17.74, 17.27,
                           16.79, 16.30, 16.01, 15.76, 15.54]
        capVolatilities = np.array(capVolatilities)/100.0
    else:
        capVolatilities = [flatVol] * numPeriods
        capVolatilities = np.array(capVolatilities)

    dayCountType = FinDayCountTypes.ACT_ACT_ISDA
    volCurve = FinCapVolCurve(valuationDate,
                              capVolDates,
                              capVolatilities,
                              dayCountType)

    zetas = np.zeros(numFwds)
    t = 0.0
    for ix in range(0, numFwds):
        t = t + dt
        zetas[ix] = volCurve.capletVol(t)

    print(zetas)
    return zetas

###############################################################################


def getForwardCurve(numFwds, r):

    fwd0 = np.zeros(numFwds)
    for i in range(0, numFwds):
        fwd0[i] = r

    fwd0 = np.array(fwd0)
    return fwd0

###############################################################################


def test_Swaptions():

    dt = 0.5
    texp = 3.0
    tmat = 10.0
    a = int(2*texp)
    b = int(2*tmat)
    numFwds = 20
    taus = np.array([dt] * numFwds)

    r = 0.05
    fwd0 = getForwardCurve(numFwds, r)
    correl = getCorrelationMatrix(numFwds, 0.00000000000001, dt)

    fwdRateVol = 0.20
    zetas = getVolCurve(numFwds, dt, fwdRateVol)

    seed = 1489
    numPaths = 2000 # 100000
    fwdsNF = LMMSimulateFwdsNF(numFwds, numPaths, fwd0,
                               zetas, correl, taus, seed)
    strike = r
    payerSwaption = 1
    useSobol = 0
    numeraireIndex = 0

    fwds1F = LMMSimulateFwds1F(numFwds, numPaths, numeraireIndex, fwd0,
                               zetas, taus, useSobol, seed)

#    LMMPrintForwards(fwds1F)   
#    LMMPrintForwards(fwdsNF)
    
    for iExp in range(1, 10):

        texp = float(iExp)
        a = int(2*texp)
        print(a, b)

        swaptionPrice1F = LMMSwaptionPricer(strike, a, b, numPaths,
                                            fwd0, fwds1F, taus, payerSwaption)

        swaptionPriceNF = LMMSwaptionPricer(strike, a, b, numPaths,
                                            fwd0, fwdsNF, taus, payerSwaption)

        swaptionVol = LMMSwaptionVolApprox(a, b, fwd0, taus, zetas, correl)

        swapVolSim1F = LMMSimSwaptionVol(a, b, fwd0, fwds1F, taus)
        swapVolSimNF = LMMSimSwaptionVol(a, b, fwd0, fwdsNF, taus)

        valuationDate = FinDate(1, 1, 2010)
        liborCurve = FinFlatCurve(valuationDate, r, 4)
        settlementDate = valuationDate
        exerciseDate = settlementDate.addMonths(a*3)
        maturityDate = settlementDate.addMonths(b*3)

        fixedCoupon = strike
        fixedFrequencyType = FinFrequencyTypes.QUARTERLY
        fixedDayCountType = FinDayCountTypes.ACT_ACT_ISDA
        floatFrequencyType = FinFrequencyTypes.QUARTERLY
        floatDayCountType = FinDayCountTypes.ACT_ACT_ISDA
        notional = 1.0

        # Pricing a PAYER
        swaptionType = FinLiborSwaptionTypes.PAYER
        swaption = FinLiborSwaption(settlementDate,
                                    exerciseDate,
                                    maturityDate,
                                    swaptionType,
                                    fixedCoupon,
                                    fixedFrequencyType,
                                    fixedDayCountType,
                                    notional,
                                    floatFrequencyType,
                                    floatDayCountType)

        model = FinModelBlack(swaptionVol)
        blackSwaptionPrice = swaption.value(valuationDate, liborCurve, model)

        print("K:%6.5f texp:%8.2f FwdVol:%9.5f SimVol1F:%9.5f SimVolNF:%9.5f RebVol:%9.5f SimPx1F:%9.5f SimPxNF:%9.5f Black Px:%9.5f" 
              % (strike, texp, fwdRateVol, swapVolSim1F, swapVolSimNF, swaptionVol, swaptionPrice1F, swaptionPriceNF, blackSwaptionPrice))

#        print(swaption)

###############################################################################


def test_CapsFloors():

    numFwds = 40
    numPaths = 100000
    dt = 0.25
    taus = np.array([dt] * numFwds)
    seeds = [42] #, 143, 101, 993, 9988]

    r = 0.04
    fwd0 = getForwardCurve(numFwds, r)

    fwdRateVol = 0.20
    zetas = getVolCurve(numFwds, dt, fwdRateVol)

    correl = getCorrelationMatrix(numFwds, 100.0, dt)

    # At the money
    K = r
    capletPricesBlack = priceCapsBlack(fwd0, zetas, numFwds, K, taus)

    numFactors = 1
    numeraireIndex = 1
    useSobol = 1

    # Examine variance for different seeds
    for seed in seeds:

        print("=============================================================")
        print("Seed:", seed)
        if numFactors == 1:
            fwds = LMMSimulateFwds1F(numFwds, numPaths, numeraireIndex, fwd0,
                                     zetas, taus, useSobol, seed)
        else:
            fwds = LMMSimulateFwdsNF(numFwds, numPaths, fwd0,
                                     zetas, correl, taus, seed)

        sumCap = LMMCapFlrPricer(numFwds, numPaths, K, fwd0, fwds, taus, 1)
        sumFlr = LMMCapFlrPricer(numFwds, numPaths, K, fwd0, fwds, taus, 0)

        print("i     CAPLET       FLRLET     BS CAP")
        for i in range(0, numFwds):
            bsCapFlrLet = capletPricesBlack[i]
            print("%d %9.6f %9.6f  %9.6f" % (i, sumCap[i]* 100.0, sumFlr[i] * 100.0,
                                             bsCapFlrLet*100.0))

###############################################################################


def test_HullBookExamples():
    ''' Examining examples on page 770 of Hull OFODS '''

    numFwds = 11
    dt = 1.00
    taus = np.array([dt] * numFwds)
    seed = 43

    r = 0.05127
    fwd0 = np.zeros(numFwds)
    for i in range(0, numFwds):
        fwd0[i] = r

    print("FWD CURVE:", fwd0)

    numPaths = 100000
    spread = 0.0025  # basis points

    ###########################################################################
    # HULL TABLE 32.1
    ###########################################################################

    useSobol = 1
    numeraireIndex = 0

    print("##################################################################")
    numFactors = 1
    lambdas1FList = [0.0, 0.1550, 0.2064, 0.1721, 0.1722, 0.1525,
                     0.1415, 0.1298, 0.1381, 0.1360, 0.1340]
    lambdas1F = np.array(lambdas1FList)

    # One factor model
    fwds1F = LMMSimulateFwds1F(numFwds, numPaths, numeraireIndex, fwd0,
                               lambdas1F, taus, useSobol, seed)

    vRatchetCaplets = LMMRatchetCapletPricer(spread, numFwds, numPaths,
                                             fwd0, fwds1F, taus) * 100.0

    hullRatchetCaplets1F = [0.0, 0.196, 0.207, 0.201, 0.194, 0.187, 0.1890, 0.172,
                            0.167, 0.160, 0.153]
    hullRatchetCaplets1F = np.array(hullRatchetCaplets1F)

#    print("Ratchet", numFactors) 
#    print(vRatchetCaplets)
#    print(hullRatchetCaplets1F)

    checkVectorDifferences(vRatchetCaplets, hullRatchetCaplets1F, 1e-2)
    
    vStickyCaplets = LMMStickyCapletPricer(spread, numFwds, numPaths,
                                           fwd0, fwds1F, taus) * 100.0

    hullStickyCaplets1F = [0.0, 0.196, 0.336, 0.412, 0.458, 0.484, 0.498, 0.502,
                           0.501, 0.497, 0.488]

#    print("Sticky", numFactors) 
#    print(vStickyCaplets)
#    print(hullStickyCaplets1F)

    checkVectorDifferences(vStickyCaplets, hullStickyCaplets1F, 1e-2)

    print("##################################################################")
    numFactors = 1
    lambdas1FList = [[0.0, 0.1550, 0.2064, 0.1721, 0.1722, 0.1525,
                      0.1415, 0.1298, 0.1381, 0.1360, 0.1340]]
    lambdas1F = np.array(lambdas1FList)

    # One factor model
    fwdsMF = LMMSimulateFwdsMF(numFwds, numFactors, numPaths, numeraireIndex,
                               fwd0, lambdas1F, taus, useSobol, seed)

    vRatchetCaplets = LMMRatchetCapletPricer(spread, numFwds, numPaths,
                                             fwd0, fwdsMF, taus) * 100.0

    hullRatchetCaplets1F = [0.0, 0.196, 0.207, 0.201, 0.194, 0.187, 0.1890, 0.172,
                            0.167, 0.160, 0.153]

#    print("Ratchet", numFactors) 
#    print(vRatchetCaplets)
#    print(hullRatchetCaplets1F)

    checkVectorDifferences(vRatchetCaplets, hullRatchetCaplets1F, 1e-2)

    vStickyCaplets = LMMStickyCapletPricer(spread, numFwds, numPaths,
                                           fwd0, fwdsMF, taus) * 100.0

    hullStickyCaplets1F = [0.0 , 0.196, 0.336, 0.412, 0.458, 0.484, 0.498, 0.502,
                           0.501, 0.497, 0.488]

    checkVectorDifferences(vStickyCaplets, hullStickyCaplets1F, 1e-2)

#    print("Sticky", numFactors) 
#    print(vStickyCaplets)
#    print(hullStickyCaplets1F)

    print("##################################################################")
    numFactors = 2
    lambdas2FList = [[0.0, 0.1410, 0.1952, 0.1678, 0.1711, 0.1525,
                      0.1406, 0.1265, 0.1306, 0.1236, 0.1163],
                     [0.0, -0.0645, -0.0670, -0.0384, -0.0196, 0.00,
                     0.0161, 0.0289, 0.0448, 0.0565, 0.0665]]
    lambdas2F = np.array(lambdas2FList)

    # Two factor model
    fwds2F = LMMSimulateFwdsMF(numFwds, numFactors, numPaths, numeraireIndex,
                               fwd0, lambdas2F, taus, useSobol, seed)

    vRatchetCaplets = LMMRatchetCapletPricer(spread, numFwds, numPaths,
                                             fwd0, fwds2F, taus) * 100.0
    hullRatchetCaplets2F = [0.0, 0.194, 0.207, 0.205, 0.198, 0.193, 0.189, 0.180,
                            0.174, 0.168, 0.162]

#    print("Ratchet", numFactors) 
#    print(vRatchetCaplets)
#    print(hullRatchetCaplets2F)

    checkVectorDifferences(vRatchetCaplets, hullRatchetCaplets2F, 1e-2)

    vStickyCaplets = LMMStickyCapletPricer(spread, numFwds, numPaths,
                                           fwd0, fwds2F, taus) * 100.0
    print(vStickyCaplets)

    hullStickyCaplets2F = [0.0, 0.196, 0.334, 0.413, 0.462, 0.492, 0.512, 0.520,
                           0.523, 0.523, 0.519]

#    print("Sticky", numFactors) 
#    print(vStickyCaplets)
#    print(hullStickyCaplets2F)

    checkVectorDifferences(vStickyCaplets, hullStickyCaplets2F, 1e-2)

    print("##################################################################")
    numFactors = 3
    lambdas3FList = [[0.0, 0.1365, 0.1928, 0.1672, 0.1698, 0.1485,
                     0.1395, 0.1261, 0.1290, 0.1197, 0.1097],
                     [0.0, -0.0662, -0.0702, -0.0406, -0.0206, 0.00,
                     0.0169, 0.0306, 0.0470, 0.0581, 0.0666],
                     [0.0, 0.0319, 0.0225, 0.000, -0.0198, -0.0347,
                     -0.0163, 0.000, 0.0151, 0.0280, 0.0384]]
    lambdas3F = np.array(lambdas3FList)

    # Three factor model
    fwds3F = LMMSimulateFwdsMF(numFwds, numFactors, numPaths, numeraireIndex,
                               fwd0, lambdas3F, taus, useSobol, seed)

    hullRatchetCaplets3F = [0.0, 0.194, 0.207, 0.205, 0.198, 0.193, 0.189, 0.180,
                            0.174, 0.168, 0.162]

    vRatchetCaplets = LMMRatchetCapletPricer(spread, numFwds, numPaths,
                                             fwd0, fwds3F, taus) * 100.0

#    print("Ratchet", numFactors) 
#    print(vRatchetCaplets)
#    print(hullRatchetCaplets3F)

    checkVectorDifferences(vRatchetCaplets, hullRatchetCaplets3F, 1e-2)

    vStickyCaplets = LMMStickyCapletPricer(spread, numFwds, numPaths,
                                           fwd0, fwds3F, taus) * 100.0

    hullStickyCaplets3F = [0.0, 0.195, 0.336, 0.418, 0.472, 0.506, 0.524, 0.533,
                           0.537, 0.537, 0.534]

#    print("Sticky", numFactors) 
#    print(vStickyCaplets)
#    print(hullStickyCaplets3F)

    checkVectorDifferences(vStickyCaplets, hullStickyCaplets3F, 1e-2)

    print("##################################################################")

###############################################################################


def test_Swap():

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

    K = r
    fwds = LMMSimulateFwdsNF(numFwds, numPaths, fwd0, zetas, correl, taus, seed)

    start = time.time()
    cpn = 0.05
    numPeriods = 10
    swap = LMMSwapPricer(cpn, numPeriods, numPaths, fwd0, fwds, taus)
    end = time.time()
    print("PRICER Period:", end - start)

    print(swap)

###############################################################################


def fwdfwdCorrelation(fwds):

    numPaths = len(fwds)
    numFwds = len(fwds[0])
    start = time.time()
    fwdCorr = LMMFwdFwdCorrelation(numFwds, numPaths, 1, fwds)
    end = time.time()
    print("CORR Period:", end - start)
    print(fwdCorr)

###############################################################################

#test_HullBookExamples()
test_CapsFloors()
# test_Swaptions()
