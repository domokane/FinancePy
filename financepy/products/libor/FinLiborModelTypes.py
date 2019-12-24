# -*- coding: utf-8 -*-

###############################################################################


class FinLiborModel(object):
    ''' This is a parent class for Libor models. '''

    def __init__(self):
        self._parentType = None
        self._volatility = 0.0
        self._implementation = 0
        pass

###############################################################################


class FinLiborModelBlack(FinLiborModel):
    def __init__(self, volatility):
        self._parentType = FinLiborModel
        self._volatility = volatility
        self._implementation = 0

###############################################################################


class FinLiborModelShiftedBlack(FinLiborModel):
    def __init__(self, volatility, shift):
        self._parentType = FinLiborModel
        self._volatility = volatility
        self._shift = shift
        self._implementation = 0

###############################################################################


class FinLiborModelSABR(FinLiborModel):
    def __init__(self, alpha, beta, rho, nu):
        self._parentType = FinLiborModel
        self._alpha = alpha
        self._beta = beta
        self._rho = rho
        self._nu = nu

###############################################################################
