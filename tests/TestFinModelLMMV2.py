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

from financepy.models.FinModelRatesLMMV2 import LMMSimulateFwds
from financepy.models.FinModelRatesLMMV2 import LMMSwaptionPricer
from financepy.models.FinModelRatesLMMV2 import LMMSwaptionVol
from financepy.models.FinModelRatesLMMV2 import LMMSwaptionVolApprox
from financepy.models.FinModelRatesLMMV2 import LMMCapFlrPricer
from financepy.models.FinModelRatesLMMV2 import priceCapsBlack, LMMSwapPricer
from financepy.models.FinModelRatesLMMV2 import LMMFwdFwdCorrelation

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
    else:
        capVolatilities = [flatVol] * numPeriods

    capVolatilities = np.array(capVolatilities)/100.0
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
    fwds = LMMSimulateFwds(numFwds, numPaths, fwd0, zetas, correl, taus, seed)

    sumSimValue = 0.0
    sumBlkValue = 0.0

    for i in range(0, len(aVector)):

        a = aVector[i]
        K = kVector[i]

        swaptionPrice = LMMSwaptionPricer(K, a, b, numPaths, fwd0, fwds, taus)
        swaptionVol = LMMSwaptionVolApprox(a, b, fwd0, taus, zetas, correl)
        swapVolSim = LMMSwaptionVol(a, b, fwd0, fwds, taus)

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
    correl = getCorrelationMatrix(numFwds, 1000.0, dt)

    fwdRateVol = 40.0
    zetas = getVolCurve(numFwds, dt, fwdRateVol)

    seed = 1489
    numPaths = 20000
    fwds = LMMSimulateFwds(numFwds, numPaths, fwd0, zetas, correl, taus, seed)
    strike = r
    payerSwaption = 1

    for iExp in range(1, 10):

        texp = float(iExp)
        a = int(4*texp)

        swaptionPrice = LMMSwaptionPricer(strike, a, b, numPaths,
                                          fwd0, fwds, taus, payerSwaption)
        swaptionVol = LMMSwaptionVolApprox(a, b, fwd0, taus, zetas, correl)
        swapVolSim = LMMSwaptionVol(a, b, fwd0, fwds, taus)

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
        v_swaption = swaption.value(valuationDate, liborCurve, model)

        print("FwdVol: %9.5f  K: %9.5f texp: %9.5f SimPx: %9.5f  SimVol:%9.5f RebVol: %9.5f Black Px: %9.5f" 
              % (fwdRateVol, strike, texp, swaptionPrice, swapVolSim,
              swaptionVol, v_swaption))

#        print(swaption)

###############################################################################


def test_CapsFloors():

    numFwds = 40
    numPaths = 10000
    dt = 0.25
    taus = np.array([dt] * numFwds)
    seeds = [42, 143, 101, 993, 9988]

    r = 0.04
    fwd0 = getForwardCurve(numFwds, r)

    fwdRateVol = 20.0
    zetas = getVolCurve(numFwds, dt, fwdRateVol)

    correl = getCorrelationMatrix(numFwds, 100.0, dt)

    # At the money
    K = r
    capletPricesBlack = priceCapsBlack(fwd0, zetas, numFwds, K, taus)

    # Examine variance for different seeds
    for seed in seeds:

        print("=============================================================")
        print("Seed:", seed)   
        fwds = LMMSimulateFwds(numFwds, numPaths, fwd0, zetas, correl, taus, seed)

        sumCap = LMMCapFlrPricer(numFwds, numPaths, K, fwd0, fwds, taus, 1)
        sumFlr = LMMCapFlrPricer(numFwds, numPaths, K, fwd0, fwds, taus, 0)

        print("i     CAP       FLR     BS CAP")
        bsCap = 0.0
        for i in range(0, numFwds):
            bsCapFlrLet = capletPricesBlack[i]
            print("%d %9.6f %9.6f  %9.6f" % (i, sumCap[i], sumFlr[i], bsCapFlrLet*100.0))

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
    fwds = LMMSimulateFwds(numFwds, numPaths, fwd0, zetas, correl, taus, seed)
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


def test_Swap():

    numFwds = 40
    numPaths = 10000
    dt = 0.25
    taus = np.array([dt] * numFwds)
    seed = 143

    r = 0.04
    fwd0 = getForwardCurve(numFwds, r)

    fwdRateVol = 40.0
    zetas = getVolCurve(numFwds, dt, fwdRateVol)
    correl = getCorrelationMatrix(numFwds, 100.0)

    K = r
    fwds = LMMSimulateFwds(numFwds, numPaths, fwd0, zetas, correl, taus, seed)

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

# test_CapsFloors()
# test_NAG_Caplet()
# test_NAG_Swaptions()
test_Swaptions()
