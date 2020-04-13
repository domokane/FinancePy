# -*- coding: utf-8 -*-


class FinFXModel(object):
    ''' This is a parent class for equity models. '''

    def __init__(self):
        self._parentType = None
        self._volatility = 0.0
        self._implementation = 0
        pass

###############################################################################


class FinFXModelBlackScholes(FinFXModel):
    def __init__(self, volatility):
        self._parentType = FinFXModel
        self._volatility = volatility
        self._implementation = 0

###############################################################################


class FinFXModelHeston(FinFXModel):
    def __init__(self, volatility, meanReversion):
        self._parentType = FinFXModel
        self._volatility = volatility
        self._meanReversion = meanReversion
        self._implementation = 0

###############################################################################
