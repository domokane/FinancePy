##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from math import sqrt
from numpy import ndarray
from numba import njit, float64, int32

##########################################################################


@njit(float64(float64[:]), fastmath=True, cache=True)
def mean(x: float):
    ''' Calculate the arithmetic mean of a vector of numbers x. '''
    n = len(x)
    m = 0.0
    for i in range(0, n):
        m += x[i]
    m = m / n
    return m

##########################################################################


@njit(float64(float64[:]), fastmath=True, cache=True)
def stdev(x: ndarray):
    ''' Calculate the standard deviation of a vector of numbers x. '''
    n = len(x)
    m = mean(x)
    v = 0.0
    for i in range(0, n):
        v += (x[i] - m) * (x[i] - m)
    sd = sqrt(v / n)
    return sd

##########################################################################


@njit(float64(float64[:]), fastmath=True, cache=True)
def stderr(x: ndarray):
    ''' Calculate the standard error estimate of a vector of numbers x. '''
    n = len(x)
    s = stdev(x)
    serr = s / sqrt(n)
    return serr

##########################################################################


@njit(float64(float64[:]), fastmath=True, cache=True)
def var(x: ndarray):
    ''' Calculate the variance of a vector of numbers x. '''
    s = stdev(x)
    v = s * s
    return v

##########################################################################


@njit(float64(float64[:], int32), fastmath=True, cache=True)
def moment(x: ndarray,
           m: int):
    ''' Calculate the m-th moment of a vector of numbers x. '''
    n = len(x)
    s = 0.0
    for i in range(0, n):
        s += pow(x[i], m)
    s = s / n
    return s

##########################################################################


@njit(float64(float64[:], float64[:]), fastmath=True, cache=True)
def correlation(x1: ndarray,
                x2: ndarray):
    ''' Calculate the correlation between two series x1 and x2. '''

    n1 = len(x1)
    n2 = len(x2)

    if n1 != n2:
        raise ValueError("Vectors have different lengths")

    m1 = mean(x1)
    m2 = mean(x2)
    sd1 = stdev(x1)
    sd2 = stdev(x2)

    prod = 0.0
    for i in range(0, n1):
        prod += x1[i] * x2[i]

    prod /= n1
    num = prod - m1 * m2
    den = sd1 * sd2
    corr = num / den
    return corr

###############################################################################
