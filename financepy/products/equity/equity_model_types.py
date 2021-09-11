##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...utils.helpers import label_to_string

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
#         s = label_to_string("OBJECT TYPE", type(self).__name__)
#         s += label_to_string("VOLATILITY", self._volatility)
#         s += label_to_string("IMPLEMENTATION", self._implementation)
#         s += label_to_string("PARAMETERS", self._parameters)
#         return s

###############################################################################


class EquityModelHeston(EquityModel):
    def __init__(self, volatility, mean_reversion):
        self._parentType = EquityModel
        self._volatility = volatility
        self._mean_reversion = mean_reversion
        self._implementation = 0

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("VOLATILITY", self._volatility)
        s += label_to_string("MEAN REVERSION", self._mean_reversion)
        s += label_to_string("IMPLEMENTATION", self._implementation)
        return s

###############################################################################
