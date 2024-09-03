##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

import sys
from typing import Union

import numpy as np
from numba import njit, float64
from prettytable import PrettyTable

from .date import Date
from .global_vars import g_days_in_year, g_small
from .error import FinError
from .day_count import DayCountTypes, DayCount


###############################################################################


def _func_name():
    """ Extract calling function name - using a protected method is not that
    advisable but calling inspect.stack is so slow it must be avoided. """
    ff = sys._getframe().f_back.f_code.co_name
    return ff

###############################################################################


def grid_index(t, grid_times):
    n = len(grid_times)
    for i in range(0, n):
        grid_time = grid_times[i]
        if abs(grid_time - t) < g_small:
            print(t, grid_times, i)
            return i

    raise FinError("Grid index not found")


###############################################################################


def beta_vector_to_corr_matrix(betas):
    """ Convert a one-factor vector of factor weights to a square correlation
    matrix. """

    num_assets = len(betas)
    correlation = np.ones(shape=(num_assets, num_assets))
    for i in range(0, num_assets):
        for j in range(0, i):
            c = betas[i] * betas[j]
            correlation[i, j] = c
            correlation[j, i] = c

    return np.array(correlation)


###############################################################################


def pv01_times(t: float,
               f: float):
    """ Calculate a bond style pv01 by calculating remaining coupon times for a
    bond with t years to maturity and a coupon frequency of f. The order of the
    list is reverse time order - it starts with the last coupon date and ends
    with the first coupon date. """

    dt = 1.0 / f
    ptimes = []

    while t >= 0.0:
        ptimes.append(t)
        t -= dt

    return ptimes


###############################################################################


def times_from_dates(dt: (Date, list),
                     value_dt: Date,
                     day_count_type: DayCountTypes = None):
    """ If a single date is passed in then return the year from valuation date
    but if a whole vector of dates is passed in then convert to a vector of
    times from the valuation date. The output is always a numpy vector of times
    which has only one element if the input is only one date. """

    if isinstance(value_dt, Date) is False:
        raise FinError("Valuation date is not a Date")

    if day_count_type is None:
        dc_counter = None
    else:
        dc_counter = DayCount(day_count_type)

    if isinstance(dt, Date):
        num_dts = 1
        times = [None]
        if dc_counter is None:
            times[0] = (dt - value_dt) / g_days_in_year
        else:
            times[0] = dc_counter.year_frac(value_dt, dt)[0]

        return times[0]

    elif isinstance(dt, list) and isinstance(dt[0], Date):
        num_dts = len(dt)
        times = []
        for i in range(0, num_dts):
            if dc_counter is None:
                t = (dt[i] - value_dt) / g_days_in_year
            else:
                t = dc_counter.year_frac(value_dt, dt[i])[0]
            times.append(t)

        return np.array(times)

    elif isinstance(dt, np.ndarray):

        raise FinError("You passed an ndarray instead of dates.")

    else:

        raise FinError("Discount factor must take dates.")

    return None

###############################################################################


def check_vector_differences(x: np.ndarray,
                             y: np.ndarray,
                             tol: float = 1e-6):
    """ Compare two vectors elementwise to see if they are more different than
    tolerance. """

    n1 = len(x)
    n2 = len(y)
    if n1 != n2:
        raise FinError("Vectors x and y do not have same size.")

    for i in range(0, n1):
        diff = x[i] - y[i]
        if abs(diff) > tol:
            print("Vector difference of:", diff, " at index: ", i)


###############################################################################


def check_dt(d: Date):
    """ Check that input d is a Date. """

    if isinstance(d, Date) is False:
        raise FinError("Should be a date dummy!")


###############################################################################


def dump(obj):
    """ Get a list of all of the attributes of a class (not built in ones) """

    attrs = dir(obj)

    non_function_attributes = [attr for attr in attrs
                               if not callable(getattr(obj, attr))]

    non_internal_attributes = [attr for attr in non_function_attributes
                               if not attr.startswith('__')]

    private_attributes = [attr for attr in non_internal_attributes
                          if attr.startswith('_')]

    public_attributes = [attr for attr in non_internal_attributes
                         if not attr.startswith('_')]

    print("PRIVATE ATTRIBUTES")
    for attr in private_attributes:
        x = getattr(obj, attr)
        print(attr, x)

    print("PUBLIC ATTRIBUTES")
    for attr in public_attributes:
        x = getattr(obj, attr)
        print(attr, x)


###############################################################################


def print_tree(array: np.ndarray,
               depth: int = None):
    """ Function that prints a binomial or trinonial tree to screen for the
    purpose of debugging. """
    n1, n2 = array.shape

    if depth is not None:
        n1 = depth

    for j in range(0, n2):
        for i in range(0, n1):
            x = array[i, n2 - j - 1]
            if x != 0.0:
                print("%10.5f" % x, end="")
            else:
                print("%10s" % '-', end="")
        print("")


###############################################################################


def input_time(dt: Date,
               curve):
    """ Validates a time input in relation to a curve. If it is a float then
    it returns a float as long as it is positive. If it is a Date then it
    converts it to a float. If it is a Numpy array then it returns the array
    as long as it is all positive. """

    small = 1e-8

    def check(t):
        if t < 0.0:
            raise FinError("Date " + str(dt) +
                           " is before curve date " + str(curve._curve_dt))
        elif t < small:
            t = small
        return t

    if isinstance(dt, float):
        t = dt
        return check(t)
    elif isinstance(dt, Date):
        t = (dt - curve.value_dt) / g_days_in_year
        return check(t)
    elif isinstance(dt, np.ndarray):
        t = dt
        if np.any(t) < 0:
            raise FinError("Date is before curve value date.")
        t = np.maximum(small, t)
        return t
    else:
        raise FinError("Unknown type.")


###############################################################################


@njit(fastmath=True, cache=True)
def listdiff(a: np.ndarray,
             b: np.ndarray):
    """ Calculate a vector of differences between two equal sized vectors. """

    if len(a) != len(b):
        raise FinError("Cannot diff lists with different sizes")

    diff = []
    for x, y in zip(a, b):
        diff.append(x - y)

    return diff


###############################################################################


@njit(fastmath=True, cache=True)
def dotproduct(x_vector: np.ndarray,
               y_vector: np.ndarray):
    """ Fast calculation of dot product using Numba. """

    dotprod = 0.0
    n = len(x_vector)
    for i in range(0, n):
        dotprod += x_vector[i] * y_vector[i]
    return dotprod

###############################################################################


@njit(fastmath=True, cache=True)
def frange(start: int,
           stop: int,
           step: int):
    """ fast range function that takes start value, stop value and step. """
    x = []
    while start <= stop:
        x.append(start)
        start += step

    return x


###############################################################################


@njit(fastmath=True, cache=True)
def normalise_weights(wt_vector: np.ndarray):
    """ Normalise a vector of weights so that they sum up to 1.0. """

    n = len(wt_vector)
    sum_wts = 0.0
    for i in range(0, n):
        sum_wts += wt_vector[i]
    for i in range(0, n):
        wt_vector[i] = wt_vector[i] / sum_wts
    return wt_vector


###############################################################################


def label_to_string(label: str,
                    value: (float, str),
                    separator: str = "\n",
                    list_format: bool = False):
    """ Format label/value pairs for a unified formatting. """
    # Format option for lists such that all values are aligned:
    # Label: value1
    #        value2
    #        ...
    label = str(label)

    if list_format and isinstance(value, list) and len(value) > 0:
        s = label + ": "
        label_spacing = " " * len(s)
        s += str(value[0])

        for v in value[1:]:
            s += "\n" + label_spacing + str(v)
        s += separator

        return s

    return f"{label}: {value}{separator}"

###############################################################################


def table_to_string(header: str,
                    value_table,
                    float_precision="10.7f"):
    """ Format a 2D array into a table-like string. """
    if (len(value_table) == 0 or type(value_table) is not list):
        print(len(value_table))
        return ""

    num_rows = len(value_table[0])

    s = header + "\n"
    for i in range(num_rows):
        for v_list in value_table:
            # isinstance is needed instead of type in case of pandas floats
            if isinstance(v_list[i], float):
                s += format(v_list[i], float_precision) + ", "
            else:
                s += str(v_list[i]) + ", "
        s = s[:-2] + "\n"

    return s[:-1]

###############################################################################


def format_table(header: (list, tuple),
                 rows: (list, tuple)):
    """ Format a 2D array into a table-like string.
    Similar to "table_to_string", but using a wrapper
    around PrettyTable to get a nice formatting. """

    t = PrettyTable(header)
    num_rows = len(header)

    if len(rows) == 0:
        print(len(rows))
        return ""

    for row in rows:
        if len(row) != num_rows:
            raise ValueError("Header and Row Size must match!")

        t.add_row(row)

    return t

###############################################################################


def to_usable_type(t):
    """ Convert a type such that it can be used with `isinstance` """
    if hasattr(t, '__origin__'):
        origin = t.__origin__
        # t comes from the `typing` module
        if origin is list:
            return (list, np.ndarray)
        elif origin is Union:
            types = t.__args__
            return tuple(to_usable_type(tp) for tp in types)
    else:
        # t is a normal type
        if t is float:
            return (int, float, np.float64)
        if isinstance(t, tuple):
            return tuple(to_usable_type(tp) for tp in t)

    return t


###############################################################################


@njit(float64(float64, float64[:], float64[:]), fastmath=True, cache=True)
def uniform_to_default_time(u, t, v):
    """ Fast mapping of a uniform random variable to a default time given a
    survival probability curve. """

    if u == 0.0:
        return 99999.0

    if u == 1.0:
        return 0.0

    num_points = len(v)
    index = 0

    for i in range(1, num_points):
        if u <= v[i - 1] and u > v[i]:
            index = i
            break

    if index == num_points + 1:
        t1 = t[num_points - 1]
        q1 = v[num_points - 1]
        t2 = t[num_points]
        q2 = v[num_points]
        lam = np.log(q1 / q2) / (t2 - t1)
        tau = t2 - np.log(u / q2) / lam
    else:
        t1 = t[index - 1]
        q1 = v[index - 1]
        t2 = t[index]
        q2 = v[index]
        tau = (t1 * np.log(q2 / u) + t2 * np.log(u / q1)) / np.log(q2 / q1)

    return tau


###############################################################################
# THIS IS NOT USED

@njit(fastmath=True, cache=True)
def accrued_tree(grid_times: np.ndarray,
                 grid_flows: np.ndarray,
                 face: float):
    """ Fast calulation of accrued interest using an Actual/Actual type of
    convention. This does not calculate according to other conventions. """

    numgrid_times = len(grid_times)

    if len(grid_flows) != numgrid_times:
        raise FinError("Grid flows not same size as grid times.")

    accrued = np.zeros(numgrid_times)

    # When the grid time is before the first coupon we have to extrapolate back

    cpn_times = np.zeros(0)
    cpn_flows = np.zeros(0)

    for i_grid in range(1, numgrid_times):

        cpn_time = grid_times[i_grid]
        cpn_flow = grid_flows[i_grid]

        if grid_flows[i_grid] > g_small:
            cpn_times = np.append(cpn_times, cpn_time)
            cpn_flows = np.append(cpn_flows, cpn_flow)

    num_cpns = len(cpn_times)

    # interpolate between coupons
    for i_grid in range(0, numgrid_times):
        t = grid_times[i_grid]
        for i in range(0, num_cpns):
            if t > cpn_times[i - 1] and t <= cpn_times[i]:
                den = cpn_times[i] - cpn_times[i - 1]
                num = t - cpn_times[i - 1]
                accrued[i_grid] = face * num * cpn_flows[i] / den
                break

    return accrued


###############################################################################


def check_argument_types(func, values):
    """ Check that all values passed into a function are of the same type
    as the function annotations. If a value has not been annotated, it
    will not be checked. """
    for value_name, annotation_type in func.__annotations__.items():

        if value_name in values:
            value = values[value_name]
            usable_type = to_usable_type(annotation_type)

        if not isinstance(value, usable_type):

            print("ERROR with function arguments for", func.__name__)
            print("This is in module", func.__module__)
            print("Please check inputs for argument >>", value_name, "<<")
            print("You have input an argument", value, "of type", type(value))
            print("The allowed types are", usable_type)
            print("It is none of these so FAILS. Please amend.")
            raise FinError("Argument Type Error")

###############################################################################
