# -*- coding: utf-8 -*-


class FinEquityModel(object):
    ''' This is a parent class for equity models. '''

    def __init__(self):
        self._parentType = None
        self._volatility = 0.0
        self._implementation = 0
        pass

###############################################################################


class FinEquityModelBlackScholes(FinEquityModel):
    def __init__(self, volatility, numStepsPerYear=100, useTree=False):
        self._parentType = FinEquityModel
        self._volatility = volatility
        self._numStepsPerYear = numStepsPerYear
        self._useTree = useTree

###############################################################################


class FinEquityModelHeston(FinEquityModel):
    def __init__(self, volatility, meanReversion):
        self._parentType = FinEquityModel
        self._volatility = volatility
        self._meanReversion = meanReversion
        self._implementation = 0

###############################################################################
