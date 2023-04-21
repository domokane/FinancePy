import numpy as np

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.frequency import annual_frequency, FrequencyTypes
from .zero_curve import BondZeroCurve
from ...market.curves.interpolator import InterpTypes
from .bond import Bond, InterpTypes, YTMCalcType


def key_rate_durations(    bond,
                           settlement_date: Date,
                           ytm: float,
                           key_rate_tenors: list = None,
                           shift: float = None,
                           rates: list = None):
            

    """
    Calculates the key rate durations for a bond.

    Parameters
    ----------
    bond : FinancePy Bond object
        
    settlement_date : FinancePy Date object
        The settlement date.
    ytm : float
        The yield to maturity.
    key_rate_tenors : list of float, optional
        The tenors of the key rates, default is None which will generate
        the tenors from 0.25 to 30 years.
    shift : float, optional
        The shift used to calculate the key rate durations, default is None
        which will set the shift to 0.0001.
    rates: list of float, optional
        Corresponding yield curve data in line with key_rate_tenors
        If None, flat yield curve is used

    Returns
    -------
    tuple of (numpy array of float, numpy array of float)
        A tuple containing the key rate tenors and the key rate durations.
    """

    # check if key_rate_tenors is None
    if key_rate_tenors is None:
    # if it is None, create an array of key rates ranging from 0.5 to 30 years
        key_rate_tenors = np.array([0.5,  1,  2,  3,  5,  7,  10, 20, 30])


    # set the shift to a small value if not give
    if not shift:
        shift = 0.0001

    # initialize an empty list for the key rate durations
    key_rate_durations = []

    # iterate over each key rate (tenor) and calculate the key rate duration
    for ind, _ in enumerate(key_rate_tenors):
        # if rates not given
        # create an array of rates where each rate is equal to the ytm value
        if rates is None:
            rates = np.ones(len(key_rate_tenors)) * ytm

        #Create set of par bonds to be used in BondZeroCurve
        # ytm and coupons are equal
        par_bonds = []
        for tenor, cpn in zip(key_rate_tenors, rates):
            mat = settlement_date.add_years(tenor)
            par_bond = Bond(settlement_date, mat, cpn,
                            bond._freq_type, bond._accrual_type)
            par_bonds.append(par_bond)

        clean_prices = [par_bond.clean_price_from_ytm(
            settlement_date, ytm, YTMCalcType.US_STREET) for par_bond, ytm in zip(par_bonds, rates)]

        par_crv = BondZeroCurve(settlement_date, par_bonds,
                                clean_prices, InterpTypes.LINEAR_ZERO_RATES)

        # calculate the full price of the bond using the discount curve
        p_zero = bond.full_price_from_discount_curve(settlement_date, par_crv)


        # shift up by the yield of corresponding par bond
        rates[ind] += shift

        par_bonds = []
        for tenor, cpn in zip(key_rate_tenors, rates):
            mat = settlement_date.add_years(tenor)
            par_bond = Bond(settlement_date, mat, cpn,
                            bond._freq_type, bond._accrual_type)
            par_bonds.append(par_bond)

        clean_prices = [par_bond.clean_price_from_ytm(
            settlement_date, ytm, YTMCalcType.US_STREET) for par_bond, ytm in zip(par_bonds, rates)]

        par_crv_up = BondZeroCurve(
            settlement_date, par_bonds, clean_prices, InterpTypes.LINEAR_ZERO_RATES)

        # calculate the full price of the bond
        # using the discount curve with the key rate shifted up
        p_up = bond.full_price_from_discount_curve(settlement_date, par_crv_up)

        # create a curve again with the key rate shifted down
        # by twice the shift value.
        rates[ind] -= shift * 2

        par_bonds = []
        for tenor, cpn in zip(key_rate_tenors, rates):
            mat = settlement_date.add_years(tenor)
            par_bond = Bond(settlement_date, mat, cpn,
                            bond._freq_type, bond._accrual_type)
            par_bonds.append(par_bond)

        clean_prices = [par_bond.clean_price_from_ytm(
            settlement_date, ytm, YTMCalcType.US_STREET) for par_bond, ytm in zip(par_bonds, rates)]
        par_crv_down = BondZeroCurve(
            settlement_date, par_bonds, clean_prices, InterpTypes.LINEAR_ZERO_RATES)

        # calculate the full price of the bond using
        p_down = bond.full_price_from_discount_curve(settlement_date, par_crv_down)

        # calculate the key rate duration
        # using the formula (P_down - P_up) / (2 * shift * P_zero)
        key_rate_duration = (p_down - p_up) / (2 * shift * p_zero)

        # append the key rate duration to the key_rate_durations list
        key_rate_durations.append(key_rate_duration)
        
    return key_rate_tenors, np.array(key_rate_durations)