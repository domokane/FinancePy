# -*- coding: utf-8 -*-

from financepy.finutils.FinDate import FinDate
from financepy.products.bonds.FinBond import FinBond
from financepy.finutils.FinFrequency import FinFrequencyTypes
from financepy.finutils.FinDayCount import FinDayCountTypes
from financepy.products.bonds.FinBondOption import FinBondOption
from financepy.products.bonds.FinBondOption import FinBondOptionTypes
from financepy.models.FinModelRatesHW import FinModelRatesHW
from financepy.models.FinModelRatesBK import FinModelRatesBK

from financepy.products.libor.FinLiborSwap import FinLiborSwap
from financepy.market.curves.FinLiborCurve import FinLiborCurve

###############################################################################


def test_FinBondCallablePuttable():

    settlementDate = FinDate(1, 12, 2019)

    ###########################################################################

    dcType = FinDayCountTypes.ACT_360
    fixedFreq = FinFrequencyTypes.SEMI_ANNUAL
    swap1 = FinLiborSwap(settlementDate, "1Y", 0.0500, fixedFreq, dcType)
    swap2 = FinLiborSwap(settlementDate, "3Y", 0.0500, fixedFreq, dcType)
    swap3 = FinLiborSwap(settlementDate, "5Y", 0.0500, fixedFreq, dcType)
    swap4 = FinLiborSwap(settlementDate, "7Y", 0.0500, fixedFreq, dcType)
    swap5 = FinLiborSwap(settlementDate, "10Y", 0.0500, fixedFreq, dcType)
    swaps = [swap1, swap2, swap3, swap4, swap5]
    discountCurve = FinLiborCurve("USD_LIBOR", settlementDate, [], [], swaps)

    ###########################################################################

    maturityDate = settlementDate.addTenor("10Y")
    coupon = 0.05
    frequencyType = FinFrequencyTypes.SEMI_ANNUAL
    accrualType = FinDayCountTypes.ACT_ACT_ICMA
    bond = FinBond(maturityDate, coupon, frequencyType, accrualType)
    expiryDate = settlementDate.addTenor("18m")
    strikePrice = 105.0
    face = 100.0

    ###########################################################################

 
 bond.calculateFlowDates(settlementDate)
    couponTimes = []
    couponFlows = []
    cpn = bond._coupon/bond._frequency
    for flowDate in bond._flowDates[1:]:
        flowTime = (flowDate - settlementDate) / gDaysInYear
        couponTimes.append(flowTime)
        couponFlows.append(cpn)
    couponTimes = np.array(couponTimes)
    couponFlows = np.array(couponFlows)

    ###########################################################################
    # Set up the call and put times and prices
    ###########################################################################

    callDates = []
    callPrices = []
    callPx = 120.0
    callDates.append(settlementDate.addTenor("5Y")); callPrices.append(callPx)
    callDates.append(settlementDate.addTenor("6Y")); callPrices.append(callPx)
    callDates.append(settlementDate.addTenor("7Y")); callPrices.append(callPx)
    callDates.append(settlementDate.addTenor("8Y")); callPrices.append(callPx)

    callTimes = []
    for dt in callDates:
        t = (dt - settlementDate) / gDaysInYear
        callTimes.append(t)

    putDates = []
    putPrices = []
    putPx = 98.0
    putDates.append(settlementDate.addTenor("5Y")); putPrices.append(putPx)
    putDates.append(settlementDate.addTenor("6Y")); putPrices.append(putPx)
    putDates.append(settlementDate.addTenor("7Y")); putPrices.append(putPx)
    putDates.append(settlementDate.addTenor("8Y")); putPrices.append(putPx)

    putTimes = []
    for dt in putDates:
        t = (dt - settlementDate) / gDaysInYear
        putTimes.append(t)

    ###########################################################################

    tmat = (maturityDate - settlementDate) / gDaysInYear
    times = np.linspace(0, tmat, 20)
    dfs = np.exp(-0.05*times)
    curve = FinDiscountCurve(settlementDate, times, dfs)

    ###########################################################################

    v1 = bond.fullPriceFromDiscountCurve(settlementDate, curve)

    sigma = 0.02  # basis point volatility
    a = 0.1

    # Test convergence
    numStepsList = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    tmat = (maturityDate - settlementDate)/gDaysInYear

    print("NUMSTEPS", "BOND_ONLY", "CALLABLE_BOND", "TIME")

    for numTimeSteps in numStepsList:
        start = time.time()
        model = FinModelRatesHullWhite(a, sigma)
        model.buildTree(tmat, numTimeSteps, times, dfs)

        v2 = model.callablePuttableBond_Tree(couponTimes, couponFlows,
                                             callTimes, callPrices,
                                             putTimes, putPrices)

        end = time.time()
        period = end-start
        print(numTimeSteps, v1, v2, period)

    if 1 == 0:
        print("RT")
        printTree(model._rt, 5)
        print("BOND")
        printTree(model._bondValues, 5)
        print("OPTION")
        printTree(model._optionValues, 5)


##############################################################################


test_FinBondCallablePuttable()
