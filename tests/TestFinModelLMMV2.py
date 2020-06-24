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

from financepy.models.FinModelRatesLMMV2 import LMMSimulateFwdsNF
from financepy.models.FinModelRatesLMMV2 import LMMSimulateFwds1F
from financepy.models.FinModelRatesLMMV2 import LMMSwaptionPricer
from financepy.models.FinModelRatesLMMV2 import LMMSimSwaptionVol
from financepy.models.FinModelRatesLMMV2 import LMMSwaptionVolApprox
from financepy.models.FinModelRatesLMMV2 import LMMCapFlrPricer
from financepy.models.FinModelRatesLMMV2 import priceCapsBlack
from financepy.models.FinModelRatesLMMV2 import LMMSwapPricer
from financepy.models.FinModelRatesLMMV2 import LMMFwdFwdCorrelation
from financepy.models.FinModelRatesLMMV2 import LMMRatchetCapletPricer
from financepy.models.FinModelRatesLMMV2 import LMMStickyCapletPricer
from financepy.models.FinModelRatesLMMV2 import printForwards

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

    return zetas

###############################################################################


def getForwardCurve(numFwds, r):

    fwd0 = np.zeros(numFwds)
    for i in range(0, numFwds):
        fwd0[i] = r

    fwd0 = np.array(fwd0)
    return fwd0

###############################################################################


def test_NAG_Swaptions():

    dt = 0.25
    b = 80
    numFwds = 80
    taus = np.array([dt] * numFwds)

    aVector = [4, 4, 4, 8, 8, 8, 20, 20, 20, 28, 28, 28, 40, 40, 40]
    kVector = [0.045, 0.05, 0.055, 0.045, 0.05, 0.055, 0.045, 0.05, 0.055,
               0.045, 0.05, 0.055, 0.045, 0.05, 0.055]

    r = 0.051
    fwd0 = getForwardCurve(numFwds, r)

    correl = getCorrelationMatrix(numFwds, 0.010, dt)

    fwdRateVol = 20.0
    zetas = getVolCurve(numFwds, dt, fwdRateVol)

    seed = 1489
    numPaths = 20000
    fwds = LMMSimulateFwdsNF(numFwds, numPaths, fwd0, zetas,
                             correl, taus, seed)

    sumSimValue = 0.0
    sumBlkValue = 0.0

    for i in range(0, len(aVector)):

        a = aVector[i]
        K = kVector[i]

        swaptionPrice = LMMSwaptionPricer(K, a, b, numPaths, fwd0, fwds, taus)
        swaptionVol = LMMSwaptionVolApprox(a, b, fwd0, taus, zetas, correl)
        swapVolSim = LMMSimSwaptionVol(a, b, fwd0, fwds, taus)

        valuationDate = FinDate(1, 1, 2010)
        liborCurve = FinFlatCurve(valuationDate, r, 4)
        settlementDate = valuationDate
        exerciseDate = settlementDate.addMonths(a*3)
        maturityDate = settlementDate.addMonths(b*3)

        fixedCoupon = K
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
        swaptionBlack = swaption.value(valuationDate, liborCurve, model)

        print("a:%d FwdVol: %9.5f K:%9.5f SimPx:%9.5f SimVol:%9.5f RebVol:%9.5f Black Px:%9.5f" 
              % (a, fwdRateVol, K, swaptionPrice, swapVolSim, swaptionVol, swaptionBlack))

        sumSimValue += swaptionPrice
        sumBlkValue += swaptionBlack

    print("Sim Sum:", sumSimValue)
    print("Sum Blk:", sumBlkValue)

###############################################################################


def test_Swaptions():

    dt = 0.25
    texp = 3.0
    tmat = 10.0
    a = int(4*texp)
    b = int(4*tmat)
    numFwds = 40
    taus = np.array([dt] * numFwds)

    r = 0.05
    fwd0 = getForwardCurve(numFwds, r)
    correl = getCorrelationMatrix(numFwds, 0.00000000000001, dt)

    fwdRateVol = 0.20
    zetas = getVolCurve(numFwds, dt, fwdRateVol)

    seed = 1489
    numPaths = 100000
    fwdsNF = LMMSimulateFwdsNF(numFwds, numPaths, fwd0,
                               zetas, correl, taus, seed)
    strike = r
    payerSwaption = 1
    fwds1F = LMMSimulateFwds1F(numFwds, numPaths, fwd0, zetas, taus, seed)

    for iExp in range(1, 10):

        texp = float(iExp)
        a = int(4*texp)

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
    seeds = [42, 143, 101, 993, 9988]

    r = 0.04
    fwd0 = getForwardCurve(numFwds, r)

    fwdRateVol = 0.20
    zetas = getVolCurve(numFwds, dt, fwdRateVol)

    correl = getCorrelationMatrix(numFwds, 100.0, dt)

    # At the money
    K = r
    capletPricesBlack = priceCapsBlack(fwd0, zetas, numFwds, K, taus)

    numFactors = 1
    # Examine variance for different seeds
    for seed in seeds:

        print("=============================================================")
        print("Seed:", seed)
        if numFactors == 1:
            fwds = LMMSimulateFwds1F(numFwds, numPaths, fwd0,
                                     zetas, taus, seed)
        else:
            fwds = LMMSimulateFwdsNF(numFwds, numPaths, fwd0,
                                     zetas, correl, taus, seed)

        sumCap = LMMCapFlrPricer(numFwds, numPaths, K, fwd0, fwds, taus, 1)
        sumFlr = LMMCapFlrPricer(numFwds, numPaths, K, fwd0, fwds, taus, 0)

        print("i     CAPLET       FLRLET     BS CAP")
        for i in range(0, numFwds):
            bsCapFlrLet = capletPricesBlack[i]
            print("%d %9.6f %9.6f  %9.6f" % (i, sumCap[i], sumFlr[i],
                                             bsCapFlrLet*100.0))

###############################################################################


def test_NAG_Caplet():

    # Example in NAG document
    numFwds = 40
    numPaths = 100000
    dt = 0.25
    taus = np.array([dt] * numFwds)
    seed = 143

    r = 0.051
    fwd0 = getForwardCurve(numFwds, r)

    fwdRateVol = 20.0
    zetas = getVolCurve(numFwds, dt, fwdRateVol)

    correl = getCorrelationMatrix(numFwds, 0.001, dt)

    K = 0.050
    capletPricesBlack = priceCapsBlack(fwd0, zetas, numFwds, K, taus)

    start = time.time()
    fwds = LMMSimulateFwdsNF(numFwds, numPaths, fwd0,
                             zetas, correl, taus, seed)
    end = time.time()
    print("SIM Period:", end - start)

    ###########################################################################

    sumCap = LMMCapFlrPricer(numFwds, numPaths, K, fwd0, fwds, taus, 1)

    print("i     CAPLET SIM (bps)   CAPLET BS (bps)")
    bsCap = 0.0
    for i in range(1, numFwds):
        bsCap = bsCap + capletPricesBlack[i]

        simCaplet = sumCap[i] - sumCap[i-1]
        print("%d %9.3f %9.3f" %
              (i, simCaplet*100, capletPricesBlack[i]*10000))

###############################################################################


def test_HullBookExamples():
    ''' Examining examples on page 770 of Hull OFODS '''

    numFwds = 11
    dt = 1.00
    taus = np.array([dt] * numFwds)
    seed = 42

    r = 0.05127
    fwd0 = np.zeros(numFwds)
    for i in range(0, numFwds):
        fwd0[i] = r

    print("FWD CURVE:", fwd0)

    fwdRateVol = None
    zetas = getVolCurve(numFwds, dt, fwdRateVol)

    print("ZETA CURVE:", zetas)

    numPaths = 1000000
    spread = 0.0025  # basis points

    # One factor model
    fwds = LMMSimulateFwds1F(numFwds, numPaths, fwd0, zetas, taus, seed)
    printForwards(fwds)

    vRatchetCaplets = LMMRatchetCapletPricer(spread, numFwds, numPaths,
                                             fwd0, fwds, taus)
    print("RATCHET CAPLETS:", vRatchetCaplets)

    vStickyCaplets = LMMStickyCapletPricer(spread, numFwds, numPaths,
                                           fwd0, fwds, taus)
    print("STICKY CAPLETS:", vStickyCaplets)

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
    fwds = LMMSimulateFwdsNF(numFwds, numPaths, fwd0, zetas,
                             correl, taus, seed)

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


test_HullBookExamples()
# test_CapsFloors()
# test_Swaptions()
# test_NAG_Caplet()
# test_NAG_Swaptions()
