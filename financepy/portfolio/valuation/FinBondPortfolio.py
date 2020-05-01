# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:51:05 2016

@author: Dominic O'Kane
"""

from ...finutils.FinDate import FinDate

from scipy import optimize

###############################################################################


class FinBondPortfolio(object):
    ''' Class for fixed coupon bonds and performing related analytics. These
    are bullet bonds which means they have regular coupon payments of a known
    size that are paid on known dates plus a payment of par at maturity.'''

    def __init__(self, settlementDate, bondList):
        ''' Create FinBondPortfolio object with a list of bond objects. '''

        self._numBonds = len(bondList)
        self._settlementDate = FinDate(1900, 1, 1)

###############################################################################
