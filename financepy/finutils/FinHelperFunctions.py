# -*- coding: utf-8 -*-
"""
Created on Sun Feb 07 14:31:53 2016

@author: Dominic O'Kane
"""

from numba import njit

##########################################################################


@njit(fastmath=True, cache=True)
def listdiff(a, b):
    ''' Calculate a vector of differences between two equal sized vectors. '''

    if len(a) != len(b):
        print("Cannot diff lists with different sizes")
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
