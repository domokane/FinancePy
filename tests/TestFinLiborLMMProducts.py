# -*- coding: utf-8 -*-

from financepy.market.volatility.FinLiborCapVolCurve import FinLiborCapVolCurve
from financepy.finutils.FinDate import FinDate
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.models.FinModelBlack import FinModelBlack
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.products.libor.FinLiborSwaption import FinLiborSwaptionTypes
from financepy.products.libor.FinLiborSwaption import FinLiborSwaption

from financepy.products.libor.FinLiborCapFloor import FinLiborCapFloorTypes

from financepy.products.libor.FinLiborLMMProducts import FinLiborLMMProducts

from financepy.products.libor.FinLiborCapFloor import FinLiborCapFloor
import numpy as np

from FinTestCases import FinTestCases, globalTestCaseMode

testCases = FinTestCases(__file__, globalTestCaseMode)


###############################################################################


def test_Swaptions():

    startYear = 2020
    endYear = 2030
    valuationDate = FinDate(1, 1, startYear)
    exerciseDate = FinDate(1, 1, 2023)
    settlementDate = valuationDate
    maturityDate = FinDate(1, 1, endYear)
    fixedCoupon = 0.04

    # DEFINE THE DISCOUNT CURVE
    discountCurve = FinDiscountCurveFlat(valuationDate,
                                 0.04,
                                 FinFrequencyTypes.ANNUAL)

    swaptionVol = 15.54

    liborSwaption = FinLiborSwaption(settlementDate,
                                     exerciseDate,
                                     maturityDate,
                                     FinLiborSwaptionTypes.PAYER,
                                     fixedCoupon,
                                     FinFrequencyTypes.ANNUAL,
                                     FinDayCountTypes.ACT_360)

    model = FinModelBlack(swaptionVol/100.0)
    v_BLK = liborSwaption.value(valuationDate, discountCurve, model)

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
        liborCurve = FinDiscountCurveFlat(valuationDate, r,
                                          FinFrequencyTypes.QUARTERLY)
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

    # Define the CAP
    # The maturity date is in 10 years so the last caplet start time is in 9
    # years which in our convention means we are modelling 10 forwards
    startYear = 2020
    endYear = 2030
    valuationDate = FinDate(1, 1, startYear)
    settlementDate = valuationDate
    capMaturityDate = FinDate(1, 1, endYear)
    frequencyType = FinFrequencyTypes.ANNUAL
    dayCountType = FinDayCountTypes.ACT_360
    capFloorRate = 0.04

    # DEFINE THE DISCOUNT CURVE
    discountCurve = FinDiscountCurveFlat(valuationDate,
                                         0.04,
                                         FinFrequencyTypes.ANNUAL)

    capVol = 15.54

    liborCap = FinLiborCapFloor(settlementDate,
                                capMaturityDate,
                                FinLiborCapFloorTypes.CAP,
                                capFloorRate,
                                None,
                                FinFrequencyTypes.ANNUAL,
                                FinDayCountTypes.ACT_360)

    model = FinModelBlack(capVol/100.0)
    v_BLK = liborCap.value(valuationDate, discountCurve, model)

    ###########################################################################
    # LMM VALUATION
    ###########################################################################

    lmmProducts = FinLiborLMMProducts(settlementDate,
                                      capMaturityDate,
                                      frequencyType,
                                      dayCountType)

    # Set up forward rate vol structure
    capVolDates = []
    capletVolTenor = "1Y"
    capletDt = valuationDate
    numForwards = endYear - startYear

    # Capvol dates has numForwards + 1 elements including today
    capVolDates.append(valuationDate)
    for i in range(0, numForwards):
        capletDt = capletDt.addTenor(capletVolTenor)
        capVolDates.append(capletDt)

    # Capvol dates has numForwards + 1 elements including zero today
    capVolatilities = [capVol] * (numForwards+1)
    capVolatilities[0] = 0.0
    capVolatilities = np.array(capVolatilities)/100.0

    dayCountType = FinDayCountTypes.ACT_ACT_ISDA
    volCurve = FinLiborCapVolCurve(valuationDate,
                                   capVolDates,
                                   capVolatilities,
                                   dayCountType)

    lambdas2FList = [[0.00, 0.1410, 0.1952, 0.1678, 0.1711, 0.1525,
                      0.1406, 0.1265, 0.1306, 0.1236],
                     [0.00, -0.0645, -0.0670, -0.0384, -0.0196, 0.00,
                     0.0161, 0.0289, 0.0448, 0.0565]]
    lambdas2F = np.array(lambdas2FList)

    # Simulate paths of future Libor rates
    numFactors = 1

    for numPaths in [10000, 20000, 50000, 100000, 200000, 400000, 1000000]:

        if numFactors == 1:
            lmmProducts.simulate1F(discountCurve, volCurve, numPaths, 0, True)
        elif numFactors == 2:
            lmmProducts.simulateMF(discountCurve, numFactors, lambdas2F,
                                   numPaths, 0, True)

        v_lmm = lmmProducts.valueCapFloor(settlementDate,
                                          capMaturityDate,
                                          FinLiborCapFloorTypes.CAP,
                                          capFloorRate,
                                          FinFrequencyTypes.ANNUAL,
                                          FinDayCountTypes.ACT_360)

        err = v_lmm - v_BLK
        print(numPaths, v_lmm, v_BLK, err)

###############################################################################

# test_CapsFloors()
# test_Swaptions()
testCases.compareTestCases()
