###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

import time
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append("..")

from financepy.finutils.FinGlobalTypes import FinOptionTypes
from financepy.products.equity.FinEquityVanillaOption import FinEquityVanillaOption
from financepy.market.curves.FinDiscountCurveFlat import FinDiscountCurveFlat
from financepy.models.FinModelBlackScholes import FinModelBlackScholes
from financepy.finutils.FinDate import FinDate

from FinTestCases import FinTestCases, globalTestCaseMode
testCases = FinTestCases(__file__, globalTestCaseMode)

###############################################################################


def test_FinNumbaNumpySpeed(useSobol):

    valueDate = FinDate(1, 1, 2015)
    expiryDate = FinDate(1, 7, 2015)
    stockPrice = 100
    volatility = 0.30
    interestRate = 0.05
    dividendYield = 0.01
    seed = 1999

    model = FinModelBlackScholes(volatility)
    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)

    useSobolInt = int(useSobol)
    
    testCases.header("NUMPATHS", "VALUE_BS", "VALUE_MC", "TIME")

    callOption = FinEquityVanillaOption(expiryDate, 100.0, 
                                        FinOptionTypes.EUROPEAN_CALL)

    value = callOption.value(valueDate, stockPrice, discountCurve,
                             dividendYield,model)

    numPoints = 20
    v_exact = [value] * numPoints

    ###########################################################################
    # DO UP TO 100K AS IT IS SLOW
    ###########################################################################

    numPathsList = np.arange(1,numPoints+1,1) * 100000

    NONUMBA_NONUMPY_v = []
    NONUMBA_NONUMPY_t = []
    
    print("PURE PYTHON")
    for numPaths in numPathsList:

        start = time.time()
        valueMC = callOption.valueMC_NONUMBA_NONUMPY(valueDate, stockPrice, discountCurve,
                                      dividendYield, model, numPaths, seed, useSobolInt)
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" % (numPaths, value, valueMC, duration))

        NONUMBA_NONUMPY_v.append(valueMC)
        NONUMBA_NONUMPY_t.append(duration+1e-10)

    NUMPY_ONLY_v = []
    NUMPY_ONLY_t = []

    print("NUMPY ONLY") 
    for numPaths in numPathsList:

        start = time.time()
        valueMC = callOption.valueMC_NUMPY_ONLY(valueDate, stockPrice, discountCurve,
                                     dividendYield, model, numPaths, seed, useSobolInt)
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" % (numPaths, value, valueMC, duration))

        NUMPY_ONLY_v.append(valueMC)
        NUMPY_ONLY_t.append(duration+1e-10)

#    speedUp = np.array(NONUMBA_NONUMPY_t)/np.array(NUMPY_ONLY_t)
#    print(NUMPY_ONLY_t)
#    print(NONUMBA_NONUMPY_t)
#    print(speedUp)

    if useSobol:
        title = "SOBOL: PURE PYTHON VS NUMPY"
    else:
        title = "PSEUDORANDOM: PURE PYTHON VS NUMPY"

    plt.figure(figsize=(8,6))
    plt.plot(numPathsList, NONUMBA_NONUMPY_t, 'o-', label="PURE PYTHON")
    plt.plot(numPathsList, NUMPY_ONLY_t, 'o-', label="NUMPY ONLY")
    plt.xlabel("Number of Paths")
    plt.ylabel("Wall Time (s)")
    plt.legend()
    plt.title(title)

    plt.figure(figsize=(8,6))
    plt.plot(numPathsList, v_exact, label="EXACT")
    plt.plot(numPathsList, NONUMBA_NONUMPY_v, 'o-', label="PURE PYTHON")
    plt.plot(numPathsList, NUMPY_ONLY_v, 'o-', label="NUMPY ONLY")
 
    plt.xlabel("Number of Paths")
    plt.ylabel("Option Value")
    plt.legend()
    plt.title(title)
   
    ###########################################################################
    # DO UP TO 10 MILLION NOW THAT WE HAVE NUMPY
    ###########################################################################
    
    numPathsList = np.arange(1,numPoints+1,1) * 1000000

    NUMPY_ONLY_v = []
    NUMPY_ONLY_t = []

    print("NUMPY ONLY") 
    for numPaths in numPathsList:

        start = time.time()
        valueMC = callOption.valueMC_NUMPY_ONLY(valueDate, stockPrice, discountCurve,
                                     dividendYield, model, numPaths, seed, useSobolInt)
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" % (numPaths, value, valueMC, duration))

        NUMPY_ONLY_v.append(valueMC)
        NUMPY_ONLY_t.append(duration)

    NUMBA_NUMPY_v = []
    NUMBA_NUMPY_t = []

    print("NUMBA+NUMPY")
    for numPaths in numPathsList:

        start = time.time()
        valueMC = callOption.valueMC_NUMPY_NUMBA(valueDate, stockPrice, discountCurve,
                                     dividendYield, model, numPaths, seed, useSobolInt)
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" % (numPaths, value, valueMC, duration))

        NUMBA_NUMPY_v.append(valueMC)
        NUMBA_NUMPY_t.append(duration)

    NUMBA_ONLY_v = []
    NUMBA_ONLY_t = []

    print("NUMBA ONLY")
    for numPaths in numPathsList:

        start = time.time()
        valueMC = callOption.valueMC_NUMBA_ONLY(valueDate, stockPrice, discountCurve,
                                     dividendYield, model, numPaths, seed, useSobolInt)
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" % (numPaths, value, valueMC, duration))

        NUMBA_ONLY_v.append(valueMC)
        NUMBA_ONLY_t.append(duration)

    NUMBA_PARALLEL_v = []
    NUMBA_PARALLEL_t = []

    print("NUMBA PARALLEL")
    for numPaths in numPathsList:

        start = time.time()
        valueMC = callOption.valueMC_NUMBA_PARALLEL(valueDate, stockPrice, discountCurve,
                                     dividendYield, model, numPaths, seed, useSobolInt)
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" % (numPaths, value, valueMC, duration))

        NUMBA_PARALLEL_v.append(valueMC)
        NUMBA_PARALLEL_t.append(duration)

#    speedUp = np.array(NUMBA_ONLY_t)/np.array(NUMBA_PARALLEL_t)
#    print("PARALLEL:", speedUp)

    ###########################################################################
    # COMPUTED USING NUMSTEPS FROM 1M to 10M
    ###########################################################################
    
    CPP_t = np.array([0.075,0.155,0.223,0.313,0.359,0.421,0.495,0.556,0.64,0.702,
                      0.765,0.841,0.923,0.982,1.05,1.125,1.195,1.261,1.333,1.408])
 
    CPP_v = np.array([9.30872,9.29576,9.29422,9.29832,9.29863,9.30153,9.2994,9.3025,9.29653,9.29875,
                      9.29897,9.29996,9.29931,9.29796,9.29784,9.2992,9.3001,9.30093,9.29876,9.29921])
 
    if useSobol:
        title = "SOBOL: COMPARING OPTIMISATIONS"
    else:
        title = "PSEUDORANDOM: COMPARING OPTIMISATIONS"

    plt.figure(figsize=(8,6))
    plt.plot(numPathsList, NUMPY_ONLY_t, 'o-', label="NUMPY ONLY")
    plt.plot(numPathsList, NUMBA_NUMPY_t, 'o-', label="NUMBA + NUMPY")
 
    plt.xlabel("Number of Paths")
    plt.ylabel("Wall Time (s)")
    plt.legend()
    plt.title(title)

    ###########################################################################

    if useSobol:
        title = "SOBOL: COMPARING OPTIMISATIONS"
    else:
        title = "PSEUDORANDOM: COMPARING OPTIMISATIONS"

    plt.figure(figsize=(8,6))
    plt.plot(numPathsList, NUMPY_ONLY_t, 'o-', label="NUMPY ONLY")
    plt.plot(numPathsList, NUMBA_NUMPY_t, 'o-', label="NUMBA + NUMPY")
    plt.plot(numPathsList, NUMBA_ONLY_t, 'o-', label="NUMBA ONLY")
    plt.plot(numPathsList, NUMBA_PARALLEL_t, 'o-', label="NUMBA PARALLEL")
 
    if useSobol == False:
        plt.plot(numPathsList, CPP_t, 'o-', label="C++")

    plt.xlabel("Number of Paths")
    plt.ylabel("Wall Time (s)")
    plt.legend()
    plt.title(title)

    ###########################################################################

    plt.figure(figsize=(8,6))
    plt.plot(numPathsList, v_exact, label="EXACT")
    plt.plot(numPathsList, NUMBA_ONLY_v, 'o-', label="NUMBA ONLY")
    plt.plot(numPathsList, CPP_v, 'o-', label="C++")

    plt.xlabel("Number of Paths")
    plt.ylabel("Option Value")
    plt.legend()
    plt.title(title)

###############################################################################

def test_FinNumbaNumbaParallel(useSobol):

    valueDate = FinDate(1, 1, 2015)
    expiryDate = FinDate(1, 7, 2015)
    stockPrice = 100
    volatility = 0.30
    interestRate = 0.05
    dividendYield = 0.01
    seed = 2021

    model = FinModelBlackScholes(volatility)
    discountCurve = FinDiscountCurveFlat(valueDate, interestRate)

    useSobolInt = int(useSobol)
    
    testCases.header("NUMPATHS", "VALUE_BS", "VALUE_MC", "TIME")

    callOption = FinEquityVanillaOption(expiryDate, 100.0, 
                                        FinOptionTypes.EUROPEAN_CALL)

    value = callOption.value(valueDate, stockPrice, discountCurve,
                             dividendYield,model)

    numPoints = 20
    v_exact = [value] * numPoints

    numPathsList = np.arange(1,numPoints+1,1) * 1000000

    NUMBA_ONLY_v = []
    NUMBA_ONLY_t = []

    print("NUMBA ONLY")
    for numPaths in numPathsList:

        start = time.time()
        valueMC = callOption.valueMC_NUMBA_ONLY(valueDate, stockPrice, discountCurve,
                                     dividendYield, model, numPaths, seed, useSobolInt)
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" % (numPaths, value, valueMC, duration))

        NUMBA_ONLY_v.append(valueMC)
        NUMBA_ONLY_t.append(duration)

    NUMBA_PARALLEL_v = []
    NUMBA_PARALLEL_t = []

    print("NUMBA PARALLEL")
    for numPaths in numPathsList:

        start = time.time()
        valueMC = callOption.valueMC_NUMBA_PARALLEL(valueDate, stockPrice, discountCurve,
                                     dividendYield, model, numPaths, seed, useSobolInt)
        end = time.time()
        duration = end - start

        print("%10d %9.5f %9.5f %9.6f" % (numPaths, value, valueMC, duration))

        NUMBA_PARALLEL_v.append(valueMC)
        NUMBA_PARALLEL_t.append(duration)

    ###########################################################################

    import matplotlib.pyplot as plt

    if useSobol:
        title = "SOBOL: NUMBA VS NUMBA + PARALLEL"
    else:
        title = "PSEUDORANDOM: NUMBA VS NUMBA + PARALLEL"

    plt.figure(figsize=(8,6))
    plt.plot(numPathsList, NUMBA_ONLY_t, 'o-', label="NUMBA ONLY")
    plt.plot(numPathsList, NUMBA_PARALLEL_t, 'o-', label="NUMBA PARALLEL")
    plt.xlabel("Number of Paths")
    plt.ylabel("Wall Time (s)")
    plt.legend()
    plt.title(title)

    plt.figure(figsize=(8,6))
    plt.plot(numPathsList, v_exact, label="EXACT")
    plt.plot(numPathsList, NUMBA_ONLY_v, 'o-', label="NUMBA ONLY")
    plt.plot(numPathsList, NUMBA_PARALLEL_v, 'o-', label="NUMBA PARALLEL")
    plt.xlabel("Number of Paths")
    plt.ylabel("Option Value")
    plt.legend()
    plt.title(title)
    
###############################################################################

if 1==0:
    test_FinNumbaNumpySpeed(False)
    test_FinNumbaNumpySpeed(True)

if 1==0:
    test_FinNumbaNumbaParallel(False)
    test_FinNumbaNumbaParallel(True)

#  testCases.compareTestCases()
