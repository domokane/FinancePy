###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

from financepy.utils.math import normcdf_slow
from financepy.utils.math import N
from financepy.utils.math import normcdf_integrate


x = -3.0


def test_normcdf1():
    result = N(x)
    assert round(result, 5) == 0.00135

def test_normcdf2():
    result = normcdf_slow(x)
    assert round(result, 5) == 0.00135

def test_normcdf_integrate():
    result = normcdf_integrate(x)
    assert round(result, 5) == 0.00135
