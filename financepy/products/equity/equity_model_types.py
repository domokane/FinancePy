##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...utils.helpers import labelToString

###############################################################################


class EquityModel:
    """ This is a parent class for equity models. """

    def __init__(self):
        pass

###############################################################################


# class EquityModelBlackScholes(EquityModel):
#     def __init__(self, 
#                  volatility: float, 
#                  implementation, parameters):

#         check_argument_types(self.__init__, locals())

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


class EquityModelHeston(EquityModel):
    def __init__(self, volatility, meanReversion):
        self._parentType = EquityModel
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
