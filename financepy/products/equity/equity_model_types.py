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

#         self.volatility = volatility
#         self.implementation = implementation
#         self.parameters = parameters

#     def __repr__(self):
#         s = label_to_string("OBJECT TYPE", type(self).__name__)
#         s += label_to_string("VOLATILITY", self.volatility)
#         s += label_to_string("IMPLEMENTATION", self.implementation)
#         s += label_to_string("PARAMETERS", self.parameters)
#         return s

###############################################################################


class EquityModelHeston(EquityModel):
    def __init__(self, volatility, mean_reversion):
        self.parentType = EquityModel
        self.volatility = volatility
        self.mean_reversion = mean_reversion
        self.implementation = 0

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("VOLATILITY", self.volatility)
        s += label_to_string("MEAN REVERSION", self.mean_reversion)
        s += label_to_string("IMPLEMENTATION", self.implementation)
        return s

###############################################################################
