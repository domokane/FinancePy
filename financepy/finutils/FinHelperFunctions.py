# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 14:31:53 2016

@author: Dominic O'Kane
"""
import numpy as np
from numba import njit
from .FinDate import FinDate
from .FinGlobalVariables import gDaysInYear

##########################################################################

def printTree(array, depth=None):
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

##########################################################################


def inputFrequency(f):
    if f in [-1, 0, 1, 2, 3, 4, 6, 12]:
        return f
    else:
        raise ValueError("Unknown frequency" + str(f))

###############################################################################


def inputTime(dt, curve):
    small = 1e-8

    def check(t):
        if t < 0.0:
            raise ValueError("Date " + str(dt) +
                             " is before curve date " + str(curve._curveDate))
        elif t < small:
            t = small
        return t

    if isinstance(dt, float):
        t = dt
        return check(t)
    elif isinstance(dt, FinDate):
        t = (dt - curve._curveDate) / gDaysInYear
        return check(t)
    elif isinstance(dt, np.ndarray):
        t = dt
        if np.any(t) < 0:
            raise ValueError("Date is before curve value date.")
        t = np.maximum(small, t)
        return t
    else:
        raise ValueError("Unknown type.")

###############################################################################

@njit(fastmath=True, cache=True)
def listdiff(a, b):
    ''' Calculate a vector of differences between two equal sized vectors. '''

    if len(a) != len(b):
        raise ValueError("Cannot diff lists with different sizes")
        return []

    diff = []
    for x, y in zip(a, b):
        diff.append(x - y)

    return diff

##########################################################################


@njit(fastmath=True, cache=True)
def dotproduct(xVector, yVector):
    dotprod = 0.0
    n = len(xVector)
    for i in range(0, n):
        dotprod += xVector[i] * yVector[i]
    return dotprod

##########################################################################


@njit(fastmath=True, cache=True)
def frange(start, stop, step):
    x = []
    while start <= stop:
        x.append(start)
        start += step

    return x

##########################################################################


@njit(fastmath=True, cache=True)
def normaliseWeights(wtVector):
    n = len(wtVector)
    sumWts = 0.0
    for i in range(0, n):
        sumWts += wtVector[i]
    for i in range(0, n):
        wtVector[i] = wtVector[i]/sumWts
    return wtVector

##########################################################################
