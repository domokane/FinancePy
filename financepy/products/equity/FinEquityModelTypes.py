##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...finutils.FinHelperFunctions import labelToString

###############################################################################


class FinEquityModel(object):
    ''' This is a parent class for equity models. '''

    def __init__(self):
        pass

###############################################################################


class FinEquityModelBlackScholes(FinEquityModel):
    def __init__(self, volatility, numStepsPerYear=100, useTree=False):
        self._parentType = FinEquityModel
        self._volatility = volatility
        self._numStepsPerYear = numStepsPerYear
        self._useTree = useTree

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("VOLATILITY", self._volatility)
        s += labelToString("NUM STEPS PER YEAR", self._numStepsPerYear)
        s += labelToString("USE TREE", self._useTree)
        return s

###############################################################################


class FinEquityModelHeston(FinEquityModel):
    def __init__(self, volatility, meanReversion):
        self._parentType = FinEquityModel
        self._volatility = volatility
        self._meanReversion = meanReversion
        self._implementation = 0

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("VOLATILITY", self._volatility)
        s += labelToString("MEAN REVERSION", self._meanReversion)
        s += labelToString("IMPLEMENTATION", self._implementation)
        return s

###############################################################################
