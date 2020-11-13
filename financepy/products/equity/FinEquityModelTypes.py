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


# class FinEquityModelBlackScholes(FinEquityModel):
#     def __init__(self, 
#                  volatility: float, 
#                  implementation, parameters):

#         checkArgumentTypes(self.__init__, locals())

#         self._volatility = volatility
#         self._implementation = implementation 
#         self._parameters = parameters

#     def __repr__(self):
#         s = labelToString("OBJECT TYPE", type(self).__name__)
#         s += labelToString("VOLATILITY", self._volatility)
#         s += labelToString("IMPLEMENTATION", self._implementation)
#         s += labelToString("PARAMETERS", self._parameters)
#         return s

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
