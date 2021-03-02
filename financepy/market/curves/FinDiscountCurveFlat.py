###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################


import numpy as np

###############################################################################

from ...utils.Date import Date
from ...utils.DayCount import FinDayCountTypes
from ...utils.Frequency import FinFrequencyTypes
from ...utils.FinHelperFunctions import labelToString
from ...utils.FinHelperFunctions import checkArgumentTypes
from ...market.curves.FinDiscountCurve import FinDiscountCurve
from ...utils.FinHelperFunctions import timesFromDates
from ...market.curves.FinInterpolator import FinInterpTypes

###############################################################################
# TODO: Do I need to add a day count to ensure rate and times are linked in
#       the correct way URGENT
###############################################################################


class FinDiscountCurveFlat(FinDiscountCurve):
    """ A very simple discount curve based on a single zero rate with its
    own specified compounding method. Hence the curve is assumed to be flat.
    It is used for quick and dirty analysis and when limited information is
    available. It inherits several methods from FinDiscountCurve. """

###############################################################################

    def __init__(self,
                 valuation_date: Date,
                 flatRate: (float, np.ndarray),
                 freq_type: FinFrequencyTypes = FinFrequencyTypes.CONTINUOUS,
                 day_count_type: FinDayCountTypes = FinDayCountTypes.ACT_ACT_ISDA):
        """ Create a discount curve which is flat. This is very useful for
        quick testing and simply requires a curve date a rate and a compound
        frequency. As we have entered a rate, a corresponding day count
        convention must be used to specify how time periods are to be measured.
        As the curve is flat, no interpolation scheme is required.
        """

        checkArgumentTypes(self.__init__, locals())

        self._valuation_date = valuation_date
        self._flatRate = flatRate
        self._freq_type = freq_type
        self._day_count_type = day_count_type

        # This is used by some inherited functions so we choose the simplest
        self._interpType = FinInterpTypes.FLAT_FWD_RATES

        # Need to set up a grid of times and discount factors
        years = np.linspace(0.0, 10.0, 41)
        dates = self._valuation_date.addYears(years)

        # Set up a grid of times and discount factors for functions
        self._dfs = self.df(dates)
        self._times = timesFromDates(dates, self._valuation_date, day_count_type)

###############################################################################

    def bump(self,
             bumpSize: float):
        """ Creates a new FinDiscountCurveFlat object with the entire curve
        bumped up by the bumpsize. All other parameters are preserved."""

        rBumped = self._flatRate + bumpSize
        discCurve = FinDiscountCurveFlat(self._valuation_date,
                                         rBumped,
                                         freq_type=self._freq_type,
                                         day_count_type=self._day_count_type)
        return discCurve

###############################################################################

    def df(self,
           dates: (Date, list)):
        """ Return discount factors given a single or vector of dates. The
        discount factor depends on the rate and this in turn depends on its
        compounding frequency and it defaults to continuous compounding. It
        also depends on the day count convention. This was set in the
        construction of the curve to be ACT_ACT_ISDA. """

        # Get day count times to use with curve day count convention
        dcTimes = timesFromDates(dates,
                                 self._valuation_date,
                                 self._day_count_type)

        dfs = self._zeroToDf(self._valuation_date,
                             self._flatRate,
                             dcTimes,
                             self._freq_type,
                             self._day_count_type)

        if isinstance(dates, Date):
            return dfs[0]
        else:
            return np.array(dfs)

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("FLAT RATE", (self._flatRate))
        s += labelToString("FREQUENCY", (self._freq_type))
        s += labelToString("DAY COUNT", (self._day_count_type))
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
