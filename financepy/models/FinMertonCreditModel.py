# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 20:00:02 2019

@author: Dominic O'Kane
"""

from ..finutils.FinMath import N
from math import sqrt, log, exp
from ...finutils.FinHelperFunctions import labelToString


def mertonCreditModelValues(assetValue,
                            bondFace,
                            timeToMaturity,
                            riskFreeRate,
                            assetGrowthRate,
                            volatility):

    lvg = assetValue / bondFace

    d1 = log(lvg) + (riskFreeRate + 0.5 *
                     volatility * volatility) * timeToMaturity
    d1 = d1 / (volatility * sqrt(timeToMaturity))
    d2 = d1 - volatility * sqrt(timeToMaturity)

    evalue = assetValue * N(d1) - bondFace * \
        exp(-riskFreeRate * timeToMaturity) * N(d2)
    dvalue = assetValue * N(-d1) + bondFace * \
        exp(-riskFreeRate * timeToMaturity) * N(d2)
    spd = -(1.0 / timeToMaturity) * log(dvalue / bondFace) - riskFreeRate

    dd = log(assetValue/bondFace)
    dd += (assetGrowthRate - (volatility**2)/2.0) * timeToMaturity
    dd = dd / volatility / sqrt(timeToMaturity)
    pd = 1.0 - N(dd)
    return (evalue, dvalue, spd, dd, pd)
