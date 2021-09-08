##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

from ...utils.date import Date
from ...utils.math import scale, test_monotonicity
from ...utils.global_vars import gDaysInYear
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.helpers import input_time
from ...utils.helpers import table_to_string
from ...market.curves.interpolator import InterpTypes, interpolate
from ...utils.error import FinError
from ...utils.frequency import annual_frequency, FrequencyTypes
from ...market.curves.discount_curve import DiscountCurve
from ...utils.helpers import label_to_string

###############################################################################


def _f(df, *args):
    curve = args[0]
    valuation_date = args[1]
    bond = args[2]
    marketCleanPrice = args[3]
    num_points = len(curve._times)
    curve._values[num_points - 1] = df
    bondDiscountPrice = bond.clean_price_from_discount_curve(
        valuation_date, curve)
    obj_fn = bondDiscountPrice - marketCleanPrice
    return obj_fn

###############################################################################


class BondZeroCurve(DiscountCurve):
    """ Class to do bootstrap exact fitting of the bond zero rate curve. """

    def __init__(self,
                 valuation_date: Date,
                 bonds: list,
                 clean_prices: list,
                 interp_type: InterpTypes = InterpTypes.FLAT_FWD_RATES):
        """ Fit a discount curve to a set of bond yields using the type of
        curve specified. """

        if len(bonds) != len(clean_prices):
            raise FinError("Num bonds does not equal number of prices.")

        self._settlement_date = valuation_date
        self._valuation_date = valuation_date
        self._bonds = bonds
        self._clean_prices = np.array(clean_prices)
        self._discount_curve = None
        self._interp_type = interp_type

        times = []
        for bond in self._bonds:
            tmat = (bond._maturity_date - self._settlement_date)/gDaysInYear
            times.append(tmat)

        times = np.array(times)
        if test_monotonicity(times) is False:
            raise FinError("Times are not sorted in increasing order")

        self._yearsToMaturity = np.array(times)

        self._bootstrap_zero_rates()

###############################################################################

    def _bootstrap_zero_rates(self):

        self._times = np.array([0.0])
        self._values = np.array([1.0])
        df = 1.0

        for i in range(0, len(self._bonds)):
            bond = self._bonds[i]
            maturity_date = bond._maturity_date
            clean_price = self._clean_prices[i]
            tmat = (maturity_date - self._settlement_date) / gDaysInYear
            argtuple = (self, self._settlement_date, bond, clean_price)
            self._times = np.append(self._times, tmat)
            self._values = np.append(self._values, df)

            optimize.newton(_f, x0=df, fprime=None, args=argtuple,
                            tol=1e-8, maxiter=100, fprime2=None)

###############################################################################

    def zero_rate(self,
                  dt: Date,
                  frequencyType: FrequencyTypes = FrequencyTypes.CONTINUOUS):
        """ Calculate the zero rate to maturity date. """
        t = input_time(dt, self)
        f = annual_frequency(frequencyType)
        df = self.df(t)

        if f == 0:  # Simple interest
            zero_rate = (1.0/df-1.0)/t
        if f == -1:  # Continuous
            zero_rate = -np.log(df) / t
        else:
            zero_rate = (df**(-1.0/t) - 1) * f
        return zero_rate

###############################################################################

    def df(self,
           dt: Date):
        t = input_time(dt, self)
        z = interpolate(t, self._times, self._values, self._interp_type.value)
        return z

###############################################################################

    def survival_prob(self,
                      dt: Date):
        t = input_time(dt, self)
        q = interpolate(t, self._times, self._values, self._interp_type.value)
        return q

###############################################################################

    def fwd(self,
            dt: Date):
        """ Calculate the continuous forward rate at the forward date. """
        t = input_time(dt, self)
        dt = 0.000001
        df1 = self.df(t)
        df2 = self.df(t+dt)
        fwd = np.log(df1/df2)/dt
        return fwd

###############################################################################

    def fwd_rate(self,
                 date1: Date,
                 date2: Date,
                 day_count_type: DayCountTypes):
        """ Calculate the forward rate according to the specified
        day count convention. """

        if date1 < self._valuation_date:
            raise FinError("Date1 before curve value date.")

        if date2 < date1:
            raise FinError("Date2 must not be before Date1")

        day_count = DayCount(day_count_type)
        year_frac = day_count.year_frac(date1, date2)[0]
        df1 = self.df(date1)
        df2 = self.df(date2)
        fwd = (df1 / df2 - 1.0) / year_frac
        return fwd

###############################################################################

    def plot(self,
             title: str):
        """ Display yield curve. """

        plt.figure(figsize=(12, 6))
        plt.title(title)
        plt.xlabel('Time to Maturity (years)')
        plt.ylabel('Zero Rate (%)')

        tmax = np.max(self._yearsToMaturity)
        t = np.linspace(0.0, int(tmax+0.5), 100)

        zero_rate = self.zero_rate(t)
        zero_rate = scale(zero_rate, 100.0)
        plt.plot(t, zero_rate, label="Zero Rate Bootstrap", marker='o')
        plt.legend(loc='lower right')
        plt.ylim((min(zero_rate)-0.3, max(zero_rate)*1.1))
        plt.grid(True)

###############################################################################

    def __repr__(self):
        # TODO
        header = "TIMES,DISCOUNT FACTORS"
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        valueTable = [self._times, self._values]
        precision = "10.7f"
        s += table_to_string(header, valueTable, precision)
        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
