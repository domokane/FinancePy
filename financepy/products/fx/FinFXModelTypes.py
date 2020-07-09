##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...finutils.FinHelperFunctions import labelToString

###############################################################################


class FinFXModel(object):

    def __init__(self):
        pass

###############################################################################


class FinFXModelBlackScholes(FinFXModel):
    ''' Basic Black Scholes model. '''

    def __init__(self, volatility):
        ''' Create Black Scholes FX model object which holds volatility. '''
        self._parentType = FinFXModel
        self._modelType = FinFXModel
        self._volatility = volatility
        self._implementation = 0

###############################################################################


class FinFXModelHeston(FinFXModel):
    ''' Heston stochastic volatility model '''

    def __init__(self, volatility, meanReversion):
        ''' Create Heston FX Model which takes in volatility and mean
        reversion. '''
        self._modelType = FinFXModel
        self._volatility = volatility
        self._meanReversion = meanReversion
        self._implementation = 0

###############################################################################


class FinFXModelSABR(FinFXModel):
    ''' SABR stochastic volatility model '''

    def __init__(self, alpha, beta, rho, nu, volatility):
        ''' Create FX Model SABR which takes alpha, beta, rho, nu and
        volatility as parameters. '''
        self._modelType = FinFXModel
        self._alpha = alpha
        self._beta = beta
        self._rho = rho
        self._nu = nu
        self._implementation = 0

###############################################################################
