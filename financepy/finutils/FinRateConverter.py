# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 10:52:29 2018

@author: Dominic O'Kane
"""
from finutils.FinError import FinError

###############################################################################


class FinRateConverter(object):
    ''' Convert rates between different compounding conventions. This is not
    used. '''

    def __init__(self, frequency):

        # we permit frequency to be entered as a string or integer
        if isinstance(frequency, int):

            if frequency == 1:
                self.name = "1Y"
                self.months = 12
            elif frequency == 2:
                self.name = "6M"
                self.months = 6
            elif frequency == 4:
                self.name = "3M"
                self.months = 3
            elif frequency == 12:
                self.name = "1M"
                self.months = 1
            else:
                raise Exception("Frequency value must be 1, 2, 4 or 12")

        elif isinstance(frequency, str):

            if frequency == "12M":
                self.months = 12
                self.name = frequency
            if frequency == "1Y":
                self.months = 12
                self.name = frequency
            elif frequency == "6M":
                self.months = 6
                self.name = frequency
            elif frequency == "3M":
                self.months = 3
                self.name = frequency
            elif frequency == "1M":
                self.months = 1
                self.name = frequency
            else:
                raise Exception("Frequency value must be 1M, 3M, 6M or 12M")

        else:
            raise FinError("Invalid frequency type.")

###############################################################################

    def __repr__(self):
        s = self.name
        return s

###############################################################################
