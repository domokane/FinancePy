##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import sys
import numpy as np
from numba import njit, float64
from typing import Union
from .FinDate import FinDate
from .FinGlobalVariables import gDaysInYear, gSmall
from .FinError import FinError
from .FinDayCount import FinDayCountTypes, FinDayCount

###############################################################################


def _funcName():
    ''' Extract calling function name - using a protected method is not that
    advisable but calling inspect.stack is so slow it must be avoided. '''
    ff = sys._getframe().f_back.f_code.co_name
    return ff

###############################################################################

def gridIndex(t, gridTimes):
    
    n = len(gridTimes)    
    for i in range(0, n):
        gridTime = gridTimes[i]
        if abs(gridTime - t) < gSmall:
            print(t, gridTimes, i) 
            return i

    raise FinError("Grid index not found")

###############################################################################


def betaVectorToCorrMatrix(betas):
    ''' Convert a one-factor vector of factor weights to a square correlation
    matrix. '''

    numAssets = len(betas)
    correlation = np.ones(shape=(numAssets, numAssets))
    for i in range(0, numAssets):
        for j in range(0, i):
            c = betas[i] * betas[j]
            correlation[i, j] = c
            correlation[j, i] = c

    return np.array(correlation)

###############################################################################


def pv01Times(t: float,
              f: float):
    ''' Calculate a bond style pv01 by calculating remaining coupon times for a
    bond with t years to maturity and a coupon frequency of f. The order of the
    list is reverse time order - it starts with the last coupon date and ends
    with the first coupon date. '''

    dt = 1.0 / f
    pv01Times = []

    while t >= 0.0:
        pv01Times.append(t)
        t -= dt

    return pv01Times

###############################################################################


def timesFromDates(dt: FinDate,
                   valuationDate: FinDate,
                   dayCountType: FinDayCountTypes = None):
    ''' If a single date is passed in then return the year from valuation date
    but if a whole vector of dates is passed in then convert to a vector of
    times from the valuation date. The output is always a numpy vector of times
    which has only one element if the input is only one date. '''

    if isinstance(valuationDate, FinDate) is False:
        raise FinError("Valuation date is not a FinDate")

    if dayCountType is None:
        dcCounter = None
    else:
        dcCounter = FinDayCount(dayCountType)

    if isinstance(dt, FinDate):
        numDates = 1
        times = [None]
        if dcCounter is None:
            times[0] = (dt - valuationDate) / gDaysInYear
        else:
            times[0] = dcCounter.yearFrac(valuationDate, dt)[0]

        return times[0]

    elif isinstance(dt, list) and isinstance(dt[0], FinDate):
        numDates = len(dt)
        times = []
        for i in range(0, numDates):
            if dcCounter is None:
                t = (dt[i] - valuationDate) / gDaysInYear
            else:
                t = dcCounter.yearFrac(valuationDate, dt[i])[0]
            times.append(t)

        return np.array(times)

    elif isinstance(dt, np.ndarray):
        raise FinError("You passed an ndarray instead of dates.")
    else:
        raise FinError("Discount factor must take dates.")

    return None

###############################################################################


def checkVectorDifferences(x: np.ndarray,
                           y: np.ndarray,
                           tol: float = 1e-6):
    ''' Compare two vectors elementwise to see if they are more different than
    tolerance. '''

    n1 = len(x)
    n2 = len(y)
    if n1 != n2:
        raise FinError("Vectors x and y do not have same size.")

    for i in range(0, n1):
        diff = x[i] - y[i]
        if abs(diff) > tol:
            print("Vector difference of:", diff, " at index: ", i)

###############################################################################


def checkDate(d: FinDate):
    ''' Check that input d is a FinDate. '''

    if isinstance(d, FinDate) is False:
        raise FinError("Should be a date dummy!")

###############################################################################


def dump(obj):
    ''' Get a list of all of the attributes of a class (not built in ones) '''

    attrs = dir(obj)

    non_function_attributes = [attr for attr in attrs
                               if not callable(getattr(obj, attr))]

    non_internal_attributes = [attr for attr in non_function_attributes
                               if not attr.startswith('__')]

    private_attributes = [attr for attr in non_internal_attributes
                          if attr.startswith('_')]

    public_attributes = [attr for attr in non_internal_attributes
                         if not attr.startswith('_')]

    print("PRIVATE ATTRIBUTES")
    for attr in private_attributes:
        x = getattr(obj, attr)
        print(attr, x)

    print("PUBLIC ATTRIBUTES")
    for attr in public_attributes:
        x = getattr(obj, attr)
        print(attr, x)

###############################################################################


def printTree(array: np.ndarray,
              depth: int = None):
    ''' Function that prints a binomial or trinonial tree to screen for the
    purpose of debugging. '''
    n1, n2 = array.shape

    if depth is not None:
        n1 = depth

    for j in range(0, n2):
        for i in range(0, n1):
            x = array[i, n2-j-1]
            if x != 0.0:
                print("%10.5f" % x, end="")
            else:
                print("%10s" % '-', end="")
        print("")

###############################################################################


def inputTime(dt: FinDate,
              curve):
    ''' Validates a time input in relation to a curve. If it is a float then
    it returns a float as long as it is positive. If it is a FinDate then it
    converts it to a float. If it is a Numpy array then it returns the array
    as long as it is all positive. '''

    small = 1e-8

    def check(t):
        if t < 0.0:
            raise FinError("Date " + str(dt) +
                           " is before curve date " + str(curve._curveDate))
        elif t < small:
            t = small
        return t

    if isinstance(dt, float):
        t = dt
        return check(t)
    elif isinstance(dt, FinDate):
        t = (dt - curve._valuationDate) / gDaysInYear
        return check(t)
    elif isinstance(dt, np.ndarray):
        t = dt
        if np.any(t) < 0:
            raise FinError("Date is before curve value date.")
        t = np.maximum(small, t)
        return t
    else:
        raise FinError("Unknown type.")

###############################################################################


@njit(fastmath=True, cache=True)
def listdiff(a: np.ndarray,
             b: np.ndarray):
    ''' Calculate a vector of differences between two equal sized vectors. '''

    if len(a) != len(b):
        raise FinError("Cannot diff lists with different sizes")

    diff = []
    for x, y in zip(a, b):
        diff.append(x - y)

    return diff

###############################################################################


@njit(fastmath=True, cache=True)
def dotproduct(xVector: np.ndarray,
               yVector: np.ndarray):
    ''' Fast calculation of dot product using Numba. '''

    dotprod = 0.0
    n = len(xVector)
    for i in range(0, n):
        dotprod += xVector[i] * yVector[i]
    return dotprod

###############################################################################


@njit(fastmath=True, cache=True)
def frange(start: int,
           stop: int,
           step: int):
    ''' fast range function that takes start value, stop value and step. '''
    x = []
    while start <= stop:
        x.append(start)
        start += step

    return x

###############################################################################


@njit(fastmath=True, cache=True)
def normaliseWeights(wtVector: np.ndarray):
    ''' Normalise a vector of weights so that they sum up to 1.0. '''

    n = len(wtVector)
    sumWts = 0.0
    for i in range(0, n):
        sumWts += wtVector[i]
    for i in range(0, n):
        wtVector[i] = wtVector[i]/sumWts
    return wtVector

###############################################################################


def labelToString(label: str,
                  value: float,
                  separator: str = "\n",
                  listFormat: bool = False):
    ''' Format label/value pairs for a unified formatting. '''
    # Format option for lists such that all values are aligned:
    # Label: value1
    #        value2
    #        ...
    label = str(label)

    if listFormat and type(value) is list and len(value) > 0:
        s = label + ": "
        labelSpacing = " " * len(s)
        s += str(value[0])

        for v in value[1:]:
            s += "\n" + labelSpacing + str(v)
        s += separator

        return s
    else:
        return f"{label}: {value}{separator}"

###############################################################################


def tableToString(header: str,
                  valueTable,
                  floatPrecision="10.7f"):
    ''' Format a 2D array into a table-like string. '''
    if (len(valueTable) == 0 or type(valueTable) is not list):
        print(len(valueTable))
        return ""

    numRows = len(valueTable[0])

    s = header + "\n"
    for i in range(numRows):
        for vList in valueTable:
            # isinstance is needed instead of type in case of pandas floats
            if (isinstance(vList[i], float)):
                s += format(vList[i], floatPrecision) + ", "
            else:
                s += str(vList[i]) + ", "
        s = s[:-2] + "\n"

    return s[:-1]

###############################################################################


def toUsableType(t):
    ''' Convert a type such that it can be used with `isinstance` '''
    if hasattr(t, '__origin__'):
        origin = t.__origin__
        # t comes from the `typing` module
        if origin is list:
            return (list, np.ndarray)
        elif origin is Union:
            types = t.__args__
            return tuple(toUsableType(tp) for tp in types)
    else:
        # t is a normal type
        if t is float:
            return (int, float, np.float64)
        if isinstance(t, tuple):
            return tuple(toUsableType(tp) for tp in t)

    return t

###############################################################################


@njit(float64(float64, float64[:], float64[:]), fastmath=True, cache=True)
def uniformToDefaultTime(u, t, v):
    ''' Fast mapping of a uniform random variable to a default time given a
    survival probability curve. '''

    if u == 0.0:
        return 99999.0

    if u == 1.0:
        return 0.0

    numPoints = len(v)
    index = 0

    for i in range(1, numPoints):
        if u <= v[i - 1] and u > v[i]:
            index = i
            break

    if index == numPoints + 1:
        t1 = t[numPoints - 1]
        q1 = v[numPoints - 1]
        t2 = t[numPoints]
        q2 = v[numPoints]
        lam = np.log(q1 / q2) / (t2 - t1)
        tau = t2 - np.log(u / q2) / lam
    else:
        t1 = t[index - 1]
        q1 = v[index - 1]
        t2 = t[index]
        q2 = v[index]
        tau = (t1 * np.log(q2 / u) + t2 * np.log(u / q1)) / np.log(q2 / q1)

    return tau

###############################################################################
# THIS IS NOT USED

@njit(fastmath=True, cache=True)
def accruedTree(gridTimes: np.ndarray,
                gridFlows: np.ndarray,
                face: float):
    ''' Fast calulation of accrued interest using an Actual/Actual type of
    convention. This does not calculate according to other conventions. '''

    numGridTimes = len(gridTimes)

    if len(gridFlows) != numGridTimes:
        raise FinError("Grid flows not same size as grid times.")

    accrued = np.zeros(numGridTimes)

    # When the grid time is before the first coupon we have to extrapolate back

    couponTimes = np.zeros(0)
    couponFlows = np.zeros(0)

    for iGrid in range(1, numGridTimes):

        cpnTime = gridTimes[iGrid]
        cpnFlow = gridFlows[iGrid]

        if gridFlows[iGrid] > gSmall:
            couponTimes = np.append(couponTimes, cpnTime)
            couponFlows = np.append(couponFlows, cpnFlow)

    numCoupons = len(couponTimes)

    # interpolate between coupons
    for iGrid in range(0, numGridTimes):
        t = gridTimes[iGrid]            
        for i in range(0, numCoupons):
            if t > couponTimes[i-1] and t <= couponTimes[i]:
                den = couponTimes[i] - couponTimes[i-1]
                num = (t - couponTimes[i-1])
                accrued[iGrid] = face * num * couponFlows[i] / den
                break
     
    return accrued

###############################################################################


def checkArgumentTypes(func, values):
    ''' Check that all values passed into a function are of the same type
    as the function annotations. If a value has not been annotated, it
    will not be checked. '''
    for valueName, annotationType in func.__annotations__.items():
        value = values[valueName]
        usableType = toUsableType(annotationType)
        if(not isinstance(value, usableType)):

            print("ERROR with function arguments for", func.__name__)
            print("This is in module", func.__module__)
            print("Please check inputs for argument >>", valueName, "<<")
            print("You have input an argument", value, "of type", type(value))
            print("The allowed types are", usableType)
            print("It is none of these so FAILS. Please amend.")
#            s = f"In {func.__module__}.{func.__name__}:\n"
#            s += f"Mismatched Types: expected a "
#            s += f"{valueName} of type '{usableType.__name__}', however"
#            s += f" a value of type '{type(value).__name__}' was given."
            raise FinError("Argument Type Error")

###############################################################################
