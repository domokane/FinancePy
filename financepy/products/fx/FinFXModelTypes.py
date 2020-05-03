# -*- coding: utf-8 -*-


class FinFXModel(object):

    def __init__(self):
        pass

###############################################################################


class FinFXModelBlackScholes(FinFXModel):
    ''' Basic Black Scholes model. '''

    def __init__(self, volatility):
        self._parentType = FinFXModel
        self._modelType = FinFXModel
        self._volatility = volatility
        self._implementation = 0

###############################################################################


class FinFXModelHeston(FinFXModel):
    ''' Heston stochastic volatility model '''

    def __init__(self, volatility, meanReversion):
        self._modelType = FinFXModel
        self._volatility = volatility
        self._meanReversion = meanReversion
        self._implementation = 0

###############################################################################


class FinFXModelSABR(FinFXModel):
    ''' SABR stochastic volatility model '''

    def __init__(self, alpha, beta, rho, nu, volatility):
        self._modelType = FinFXModel
        self._alpha = alpha
        self._beta = beta
        self._rho = rho
        self._nu = nu
        self._implementation = 0

###############################################################################
