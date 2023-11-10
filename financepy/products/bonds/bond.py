###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

###############################################################################
# TODO: - ROUNDING CONVENTIONS FOR ACCRUED
# TODO: - CHECK OAS CALCULATION
# TODO:  - Check how first coupon on floating leg is sized on asset swaps. """
###############################################################################

# https://www.dmo.gov.uk/media/15004/convention_changes.pdf

###############################################################################
# Conventions:
#  GILTS - SEMI ANNUAL ACT/ACT
#  US TREASURIES
###############################################################################

###############################################################################
# NOTE THAT I ASSUME THAT IF YOU SETTLE A SWAP ON A COUPON PAYMENT DATE YOU
# GET THE COUPON AND THE ACCRUED INTEREST EQUALS THE COUPON.
###############################################################################

import numpy as np
from scipy import optimize

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.frequency import annual_frequency, FrequencyTypes
from ...utils.global_vars import gDaysInYear, gSmall
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.schedule import Schedule
from ...utils.calendar import Calendar
from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes
from ...utils.calendar import DateGenRuleTypes
from ...utils.helpers import label_to_string, check_argument_types
from ...utils.math import npv
from ...market.curves.discount_curve import DiscountCurve
from ...market.curves.interpolator import InterpTypes
from .zero_curve import BondZeroCurve


# References https://www.dmo.gov.uk/media/15011/yldeqns_v1.pdf
# DO TRUE YIELD
# JAPANESE SIMPLE YIELD

###############################################################################


from enum import Enum


class YTMCalcType(Enum):
    ZERO = 0,
    UK_DMO = 1,
    US_STREET = 2,
    US_TREASURY = 3,
    CFETS = 4  # China Foreign Exchange Trade System


###############################################################################


def _f(ytm, *args):
    """ Function used to do root search in price to yield calculation. """
    bond = args[0]
    settlement_date = args[1]
    price = args[2]
    convention = args[3]
    px = bond.dirty_price_from_ytm(settlement_date, ytm, convention)
    obj_fn = px - price
    return obj_fn

###############################################################################


def _g(oas, *args):
    """ Function used to do root search in price to OAS calculation. """
    bond = args[0]
    settlement_date = args[1]
    price = args[2]
    discount_curve = args[3]
    px = bond.dirty_price_from_oas(settlement_date, discount_curve, oas)
    obj_fn = px - price
    return obj_fn

###############################################################################


class Bond:
    """ Class for fixed coupon bonds and performing related analytics. These
    are bullet bonds which means they have regular coupon payments of a known
    size that are paid on known dates plus a payment of par at maturity. """

    def __init__(self,
                 issue_date: Date,
                 maturity_date: Date,
                 coupon: float,  # Annualised bond coupon
                 freq_type: FrequencyTypes,
                 accrual_type: DayCountTypes,
                 ex_div_days: int = 0,
                 calendar_type: CalendarTypes=CalendarTypes.WEEKEND,
                 bus_day_rule_type=BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type=DateGenRuleTypes.BACKWARD):
        """ Create Bond object by providing the issue date, maturity Date,
        coupon frequency, annualised coupon, the accrual convention type, face
        amount and the number of ex-dividend days. A calendar type is used
        to determine holidays from which coupon dates might be shifted."""

        check_argument_types(self.__init__, locals())

        if issue_date >= maturity_date:
            raise FinError("Issue Date must preceded maturity date.")

        # If the maturity date falls on the last day of the month we assume
        # that earlier flows also fall on month ends
        self._end_of_month = False
        if maturity_date.is_eom():
            self._end_of_month = True

        self._issue_date = issue_date
        self._maturity_date = maturity_date

        if coupon == 0.0 or freq_type == FrequencyTypes.ZERO or accrual_type == accrual_type.ZERO:
            raise FinError(f"Zero coupon bonds must use BondZero class. coupon:{coupon}. freq:{freq_type}" +
                           f"accrual_type:{accrual_type}")

        self._coupon = coupon
        self._freq_type = freq_type

        self._accrual_type = accrual_type
        self._frequency = annual_frequency(freq_type)

        if ex_div_days > 90:
            raise FinError("Ex dividend days cannot be more than 90"
                           + str(ex_div_days))

        self._ex_div_days = ex_div_days

        self._par = 100.0  # This is how price is quoted and amount at maturity
        self._calendar_type = calendar_type

        self._coupon_dates = []  # can be holidays or weekend
        self._payment_dates = []  # Actual pay dates are adjusted to bus days
        self._flow_amounts = []

        self._accrued_interest = None
        self._accrued_days = 0.0
        self._alpha = 0.0

        self.bus_day_rule_type = bus_day_rule_type
        self.date_gen_rule_type = date_gen_rule_type

        self._calculate_coupon_dates()
        self._calculate_flows()

    ###########################################################################

    def _calculate_coupon_dates(self):
        """ Determine the bond coupon dates. Note that for analytical
        calculations these are not usually adjusted and so may fall on a
        weekend or holiday.
        """

        # This should only be called once from init
        # bus_day_rule_type = BusDayAdjustTypes.FOLLOWING
        # date_gen_rule_type = DateGenRuleTypes.BACKWARD

        self._coupon_dates = Schedule(self._issue_date,
                                      self._maturity_date,
                                      self._freq_type,
                                      CalendarTypes.NONE,
                                      self.bus_day_rule_type,
                                      self.date_gen_rule_type,
                                      end_of_month=self._end_of_month)._generate()

    ###########################################################################

    def _calculate_payment_dates(self):
        """ For the actual payment dates, they are adjusted
        and so we then use the calendar payment dates. Although payments are
        calculated as though coupon periods are the same length, payments that
        fall on a Saturday or Sunday can only be made on the next business day
        """

        # dates are adjusted forward to the next business day
        bus_day_adj_type = BusDayAdjustTypes.FOLLOWING
        calendar = Calendar(self._calendar_type)

        self._calculate_coupon_dates()

        self._payment_dates = []

        # Expect at least an issue date and a maturity date - if not - problem
        if len(self._coupon_dates) < 2:
            raise FinError("Unable to calculate pmnt dates with only one pmnt")

        # I do not adjust the first date as it is the issue date
        self._payment_dates.append(self._coupon_dates[0])

        for cpn_date in self._coupon_dates[1:]:
            pmt_date = calendar.adjust(cpn_date,
                                       bus_day_adj_type)

            self._payment_dates.append(pmt_date)

    ###########################################################################

    def _calculate_flows(self):
        """ Determine the bond cash flow payment amounts without principal.
        There is no adjustment based on the adjusted payment dates. """

        self._flow_amounts = [0.0]

        for _ in self._coupon_dates[1:]:
            cpn = self._coupon / self._frequency
            self._flow_amounts.append(cpn)

    ###########################################################################

    def dirty_price_from_ytm(self,
                             settlement_date: Date,
                             ytm: float,
                             convention: YTMCalcType = YTMCalcType.UK_DMO):
        """ Calculate the dirty price of bond from its yield to maturity. This
        function is vectorised with respect to the yield input. It implements
        a number of standard conventions for calculating the YTM. """

        if convention not in YTMCalcType:
            raise FinError("Yield convention unknown." + str(convention))

        # TODO check that no unneccesary calculations are being done
        self.accrued_interest(settlement_date, 1.0)

        #######################################################################
        # HANDLE EX_DIVIDEND DATES
        #######################################################################

        pay_first_coupon = 1.0
        if settlement_date > self._ex_div_date:
            pay_first_coupon = 0.0

        #######################################################################

        ytm = np.array(ytm)  # VECTORIZED
        ytm = ytm + 0.000000000012345  # SNEAKY LOW-COST TRICK TO AVOID y=0

        f = annual_frequency(self._freq_type)
        c = self._coupon

        if convention == YTMCalcType.ZERO:
            raise FinError("Zero coupon bonds must use BondZero class.")

        v = 1.0 / (1.0 + ytm / f)

        # n is the number of flows after the next coupon
        n = 0
        for dt in self._coupon_dates:
            if dt > settlement_date:
                n += 1
        n = n - 1

        if n < 0:
            raise FinError("No coupons left")

        if convention == YTMCalcType.UK_DMO:
            if n == 0:
                dp = (v ** (self._alpha)) * (1.0 + pay_first_coupon * c / f)
            else:
                term1 = (c / f) * pay_first_coupon
                term2 = (c / f) * v
                term3 = (c / f) * v * v * (1.0 - v ** (n - 1)) / (1.0 - v)
                term4 = (v ** n)
                dp = (v ** (self._alpha)) * (term1 + term2 + term3 + term4)
#                print(term1, term2, term3, term4, v, self._alpha, dp)
        elif convention == YTMCalcType.US_TREASURY:
            if n == 0:
                dp = (v ** (self._alpha)) * (1.0 + c / f)
            else:
                term1 = (c / f) * pay_first_coupon
                term2 = (c / f) * v
                term3 = (c / f) * v * v * (1.0 - v ** (n - 1)) / (1.0 - v)
                term4 = (v ** n)
                vw = 1.0 / (1.0 + self._alpha * ytm / f)
                dp = (vw) * (term1 + term2 + term3 + term4)
        elif convention == YTMCalcType.US_STREET:
            if n == 0:
                vw = 1.0 / (1.0 + self._alpha * ytm / f)
                dp = vw * (1.0 + c / f)
            else:
                term1 = (c / f) * pay_first_coupon
                term2 = (c / f) * v
                term3 = (c / f) * v * v * (1.0 - v ** (n - 1)) / (1.0 - v)
                term4 = (v ** n)
                dp = (v ** (self._alpha)) * (term1 + term2 + term3 + term4)
        elif convention == YTMCalcType.CFETS:
            if n == 0:
                last_year = self._maturity_date.add_tenor("-12M")

                dc = DayCount(DayCountTypes.ACT_365L)

                alpha = (1 - dc.year_frac(last_year,
                                          settlement_date,
                                          self._maturity_date,
                                          freq_type=FrequencyTypes.ANNUAL)[0])

                vw = 1.0 / (1.0 + alpha * ytm)
                dp = vw * (1.0 + c / f)
            else:
                term1 = (c / f) * pay_first_coupon
                term2 = (c / f) * v
                term3 = (c / f) * v * v * (1.0 - v ** (n - 1)) / (1.0 - v)
                term4 = (v ** n)
                dp = (v ** (self._alpha)) * (term1 + term2 + term3 + term4)
        else:
            raise FinError("Unknown yield convention")

        return dp * self._par

    ###########################################################################

    def principal(self,
                  settlement_date: Date,
                  ytm: float,
                  face: (float),
                  convention: YTMCalcType):
        """ Calculate the principal value of the bond based on the face
        amount from its discount margin and making assumptions about the
        future Ibor rates. """

        dirty_price = self.dirty_price_from_ytm(settlement_date,
                                                ytm,
                                                convention)

        principal = dirty_price * face / self._par
        principal = principal - self._accrued_interest
        return principal

    ###########################################################################

    def dollar_duration(self,
                        settlement_date: Date,
                        ytm: float,
                        convention: YTMCalcType = YTMCalcType.UK_DMO):
        """ Calculate the risk or dP/dy of the bond by bumping. This is also
        known as the DV01 in Bloomberg. """

        dy = 0.0001  # 1 basis point
        p0 = self.dirty_price_from_ytm(settlement_date, ytm - dy, convention)
        p2 = self.dirty_price_from_ytm(settlement_date, ytm + dy, convention)
        durn = -(p2 - p0) / dy / 2.0
        return durn

    ###########################################################################

    def macauley_duration(self,
                          settlement_date: Date,
                          ytm: float,
                          convention: YTMCalcType = YTMCalcType.UK_DMO):
        """ Calculate the Macauley duration of the bond on a settlement date
        given its yield to maturity. """

        dd = self.dollar_duration(settlement_date, ytm, convention)
        fp = self.dirty_price_from_ytm(settlement_date, ytm, convention)
        md = dd * (1.0 + ytm / self._frequency) / fp
        return md

    ###########################################################################

    def modified_duration(self,
                          settlement_date: Date,
                          ytm: float,
                          convention: YTMCalcType = YTMCalcType.UK_DMO):
        """ Calculate the modified duration of the bond on a settlement date
        given its yield to maturity. """

        dd = self.dollar_duration(settlement_date, ytm, convention)
        fp = self.dirty_price_from_ytm(settlement_date, ytm, convention)
        md = dd / fp
        return md

    ###########################################################################

    def key_rate_durations(bond,
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
        # if it is None, create an array of key rates from 0.5 to 30 years

        if key_rate_tenors is None:
            key_rate_tenors = np.array([0.5,  1,  2,  3,  5,  7,  10, 20, 30])

        # set the shift to a small value if not give
        if not shift:
            shift = 0.0001

        lin_zero_interp = InterpTypes.LINEAR_ZERO_RATES
        us_street = YTMCalcType.US_STREET

        # initialize an empty list for the key rate durations
        key_rate_durations = []

        # iterate over each key rate (tenor) and calculate key rate duration
        for ind, _ in enumerate(key_rate_tenors):

            # if rates are not given create an array of rates where each
            # rate is equal to the ytm value
            if rates is None:
                rates = np.ones(len(key_rate_tenors)) * ytm

            # Create set of par bonds to be used in BondZeroCurve
            # ytm and coupons are equal
            par_bonds = []

            for tenor, cpn in zip(key_rate_tenors, rates):

                mat_date = settlement_date.add_years(tenor)

                par_bond = Bond(settlement_date, mat_date, cpn,
                                bond._freq_type, bond._accrual_type)

                par_bonds.append(par_bond)

            clean_prices = []

            for par_bond, ytm in zip(par_bonds, rates):

                clean_price = par_bond.clean_price_from_ytm(settlement_date,
                                                            ytm,
                                                            us_street)
                clean_prices.append(clean_price)

            par_crv = BondZeroCurve(settlement_date,
                                    par_bonds,
                                    clean_prices,
                                    lin_zero_interp)

            # calculate the dirty price of the bond using the discount curve
            p_zero = bond.dirty_price_from_discount_curve(settlement_date, par_crv)

            # shift up by the yield of corresponding par bond
            rates[ind] += shift

            par_bonds = []

            for tenor, cpn in zip(key_rate_tenors, rates):

                mat = settlement_date.add_years(tenor)

                par_bond = Bond(settlement_date, mat, cpn,
                                bond._freq_type, bond._accrual_type)

                par_bonds.append(par_bond)

            clean_prices = []

            for par_bond, ytm in zip(par_bonds, rates):

                clean_price = par_bond.clean_price_from_ytm(settlement_date,
                                                            ytm,
                                                            us_street)
                clean_prices.append(clean_price)

            par_crv_up = BondZeroCurve(settlement_date,
                                       par_bonds,
                                       clean_prices,
                                       lin_zero_interp)

            # calculate the full price of the bond
            # using the discount curve with the key rate shifted up
            p_up = bond.dirty_price_from_discount_curve(settlement_date, par_crv_up)

            # create a curve again with the key rate shifted down
            # by twice the shift value.
            rates[ind] -= shift * 2

            par_bonds = []

            for tenor, cpn in zip(key_rate_tenors, rates):

                mat = settlement_date.add_years(tenor)

                par_bond = Bond(settlement_date, mat, cpn,
                                bond._freq_type, bond._accrual_type)

                par_bonds.append(par_bond)

            clean_prices = []

            for par_bond, ytm in zip(par_bonds, rates):

                clean_price = par_bond.clean_price_from_ytm(settlement_date,
                                                            ytm,
                                                            us_street)
                clean_prices.append(clean_price)

            par_crv_down = BondZeroCurve(settlement_date,
                                         par_bonds,
                                         clean_prices,
                                         lin_zero_interp)

            # calculate the full price of the bond using
            p_down = bond.dirty_price_from_discount_curve(settlement_date, par_crv_down)

            # calculate the key rate duration
            # using the formula (P_down - P_up) / (2 * shift * P_zero)
            key_rate_duration = (p_down - p_up) / (2 * shift * p_zero)

            # append the key rate duration to the key_rate_durations list
            key_rate_durations.append(key_rate_duration)

        return key_rate_tenors, np.array(key_rate_durations)

    ###########################################################################

    def convexity_from_ytm(self,
                           settlement_date: Date,
                           ytm: float,
                           convention: YTMCalcType = YTMCalcType.UK_DMO):
        """ Calculate the bond convexity from the yield to maturity. This
        function is vectorised with respect to the yield input. """

        dy = 0.0001  # 1 basis point
        p0 = self.dirty_price_from_ytm(settlement_date, ytm - dy, convention)
        p1 = self.dirty_price_from_ytm(settlement_date, ytm, convention)
        p2 = self.dirty_price_from_ytm(settlement_date, ytm + dy, convention)
        conv = ((p2 + p0) - 2.0 * p1) / dy / dy / p1 / self._par
        return conv

    ###########################################################################

    def clean_price_from_ytm(self,
                             settlement_date: Date,
                             ytm: float,
                             convention: YTMCalcType = YTMCalcType.UK_DMO):
        """ Calculate the bond clean price from the yield to maturity. This
        function is vectorised with respect to the yield input. """

        dp = self.dirty_price_from_ytm(settlement_date, ytm, convention)
        accrued = self.accrued_interest(settlement_date, self._par)
        cp = dp - accrued
        return cp

    ###########################################################################

    def clean_price_from_discount_curve(self,
                                        settlement_date: Date,
                                        discount_curve: DiscountCurve):
        """ Calculate the clean bond value using some discount curve to
        present-value the bond's cash flows back to the curve anchor date and
        not to the settlement date. """

        self.accrued_interest(settlement_date, self._par)
        dirty_price = self.dirty_price_from_discount_curve(settlement_date,
                                                           discount_curve)

        accrued = self.accrued_interest(settlement_date, self._par)
        clean_price = dirty_price - accrued
        return clean_price

    ###########################################################################

    def dirty_price_from_discount_curve(self,
                                        settlement_date: Date,
                                        discount_curve: DiscountCurve):
        """ Calculate the bond price using a provided discount curve to PV the
        bond's cash flows to the settlement date. As such it is effectively a
        forward bond price if the settlement date is after the valuation date.
        """

        if settlement_date < discount_curve._valuation_date:
            raise FinError("Bond settles before Discount curve date")

        if settlement_date > self._maturity_date:
            raise FinError("Bond settles after it matures.")

        self._calc_pcd_ncd(settlement_date)

        cal = Calendar(self._calendar_type)

        self._ex_div_date = cal.add_business_days(self._ncd,
                                                  -self._ex_div_days)

        pay_first_coupon = 1.0
        if settlement_date > self._ex_div_date:
            pay_first_coupon = 0.0

        px = 0.0
        df = 1.0
        dfSettle = discount_curve.df(settlement_date)

        dt = self._coupon_dates[1]
        if dt > settlement_date:
            df = discount_curve.df(dt)
            flow = self._coupon / self._frequency
            pv = flow * df
            px += pv * pay_first_coupon

        for dt in self._coupon_dates[2:]:

            # coupons paid on a settlement date are paid to the seller
            if dt > settlement_date:
                df = discount_curve.df(dt)
                flow = self._coupon / self._frequency
                pv = flow * df
                px += pv

        px += df
        px = px / dfSettle

        return px * self._par

    ###########################################################################

    def current_yield(self, clean_price):
        """ Calculate the current yield of the bond which is the
        coupon divided by the clean price (not the full price)"""

        y = self._coupon * self._par / clean_price
        return y

    ###########################################################################

    def yield_to_maturity(self,
                          settlement_date: Date,
                          clean_price: float,
                          convention: YTMCalcType = YTMCalcType.US_TREASURY):
        """ Calculate the bond's yield to maturity by solving the price
        yield relationship using a one-dimensional root solver. """

        if type(clean_price) is float or type(clean_price) is np.float64:
            clean_prices = np.array([clean_price])
        elif type(clean_price) is list or type(clean_price) is np.ndarray:
            clean_prices = np.array(clean_price)
        else:
            raise FinError("Unknown type for clean_price "
                           + str(type(clean_price)))

        self.accrued_interest(settlement_date, 1.0)

        accrued_amount = self._accrued_interest * self._par
        dirty_prices = (clean_prices + accrued_amount)
        ytms = []

        for dirty_price in dirty_prices:

            argtuple = (self, settlement_date, dirty_price, convention)

            ytm = optimize.newton(_f,
                                  x0=0.05,  # guess initial value of 5%
                                  fprime=None,
                                  args=argtuple,
                                  tol=1e-8,
                                  maxiter=50,
                                  fprime2=None)

            ytms.append(ytm)

        if len(ytms) == 1:
            return ytms[0]
        else:
            return np.array(ytms)

    ###########################################################################

    def _calc_pcd_ncd(self,
                      settlement_date: Date):

        num_flows = len(self._coupon_dates)

        if num_flows == 0:
            raise FinError("Accrued interest - not enough flow dates.")

        for iFlow in range(1, num_flows):
            # coupons paid on a settlement date are paid to the seller
            if self._coupon_dates[iFlow] > settlement_date:
                self._pcd = self._coupon_dates[iFlow - 1]
                self._ncd = self._coupon_dates[iFlow]
                break

    ###########################################################################

    def accrued_interest(self,
                         settlement_date: Date,
                         face: float = 100.0):
        ''' Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date. Note that for some day
        count schemes (such as 30E/360) this is not actually the number of days
        between the previous coupon payment date and settlement date. If the
        bond trades with ex-coupon dates then you need to use the number of
        days before the coupon date the ex-coupon date is. You can specify the
        calendar to be used in the bond constructor - NONE means only calendar
        days, WEEKEND is only weekends or you can specify a country calendar
        for business days.'''

        self._calc_pcd_ncd(settlement_date)

        dc = DayCount(self._accrual_type)
        cal = Calendar(self._calendar_type)

        self._ex_div_date = cal.add_business_days(self._ncd,
                                                  -self._ex_div_days)

        (acc_factor, num, _) = dc.year_frac(self._pcd,
                                            settlement_date,
                                            self._ncd,
                                            self._freq_type)

        self._alpha = 1.0 - acc_factor * self._frequency

        if settlement_date > self._ex_div_date:
            self._accrued_interest = acc_factor - 1.0 / self._frequency
        else:
            self._accrued_interest = acc_factor

        self._accrued_interest *= self._coupon * face

        self._accrued_days = num

        return self._accrued_interest

    ###########################################################################

    def asset_swap_spread(
            self,
            settlement_date: Date,
            clean_price: float,
            discount_curve: DiscountCurve,
            swapFloatDayCountConventionType=DayCountTypes.ACT_360,
            swapFloatFrequencyType=FrequencyTypes.SEMI_ANNUAL,
            swapFloatCalendarType=CalendarTypes.WEEKEND,
            swapFloatBusDayAdjustRuleType=BusDayAdjustTypes.FOLLOWING,
            swapFloatDateGenRuleType=DateGenRuleTypes.BACKWARD):
        """ Calculate the par asset swap spread of the bond. The discount curve
        is a Ibor curve that is passed in. This function is vectorised with
        respect to the clean price. """

        clean_price = np.array(clean_price)
        self.accrued_interest(settlement_date, 1.0)
        accrued_amount = self._accrued_interest * self._par
        bondPrice = clean_price + accrued_amount
        # Calculate the price of the bond discounted on the Ibor curve
        pvIbor = 0.0
        prev_date = self._pcd

        for dt in self._coupon_dates[1:]:

            # coupons paid on a settlement date are paid to the seller
            if dt > settlement_date:
                df = discount_curve.df(dt)
                pvIbor += df * self._coupon / self._frequency

        pvIbor += df

        # Calculate the PV01 of the floating leg of the asset swap
        # I assume here that the coupon starts accruing on the settlement date
        prev_date = self._pcd
        schedule = Schedule(settlement_date,
                            self._maturity_date,
                            swapFloatFrequencyType,
                            swapFloatCalendarType,
                            swapFloatBusDayAdjustRuleType,
                            swapFloatDateGenRuleType)

        day_count = DayCount(swapFloatDayCountConventionType)

        prev_date = self._pcd
        pv01 = 0.0
        for dt in schedule._adjusted_dates[1:]:
            df = discount_curve.df(dt)
            year_frac = day_count.year_frac(prev_date, dt)[0]
            pv01 = pv01 + year_frac * df
            prev_date = dt

        asw = (pvIbor - bondPrice / self._par) / pv01
        return asw

    ###########################################################################

    def dirty_price_from_oas(self,
                             settlement_date: Date,
                             discount_curve: DiscountCurve,
                             oas: float):
        """ Calculate the full price of the bond from its OAS given the bond
        settlement date, a discount curve and the oas as a number. """

        self.accrued_interest(settlement_date, 1.0)
        f = self._frequency
        c = self._coupon

        pv = 0.0
        for dt in self._coupon_dates[1:]:

            # coupons paid on a settlement date are paid to the seller
            if dt > settlement_date:
                t = (dt - settlement_date) / gDaysInYear

                t = np.maximum(t, gSmall)

                df = discount_curve.df(dt)

                # determine the Ibor implied zero rate
                r = f * (np.power(df, -1.0 / t / f) - 1.0)

                # determine the OAS adjusted zero rate
                df_adjusted = np.power(1.0 + (r + oas) / f, -t * f)
                pv = pv + (c / f) * df_adjusted

        pv = pv + df_adjusted
        pv *= self._par
        return pv

    ###########################################################################

    def option_adjusted_spread(self,
                               settlement_date: Date,
                               clean_price: float,
                               discount_curve: DiscountCurve):
        """ Return OAS for bullet bond given settlement date, clean bond price
        and the discount relative to which the spread is to be computed. """

        if type(clean_price) is float or type(clean_price) is np.float64:
            clean_prices = np.array([clean_price])
        elif type(clean_price) is list or type(clean_price) is np.ndarray:
            clean_prices = np.array(clean_price)
        else:
            raise FinError("Unknown type for clean_price "
                           + str(type(clean_price)))

        self.accrued_interest(settlement_date, 1.0)

        accrued_amount = self._accrued_interest * self._par
        dirty_prices = clean_prices + accrued_amount

        oass = []

        for dirty_price in dirty_prices:
            argtuple = (self, settlement_date, dirty_price, discount_curve)

            oas = optimize.newton(_g,
                                  x0=0.01,  # initial value of 1%
                                  fprime=None,
                                  args=argtuple,
                                  tol=1e-8,
                                  maxiter=50,
                                  fprime2=None)

            oass.append(oas)

        if len(oass) == 1:
            return oass[0]
        else:
            return np.array(oass)

    ###########################################################################

    def bond_payments(self, settlement_date: Date, face: (float)):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """

        flow = face * self._coupon / self._frequency

        flow_str = ""

        for dt in self._coupon_dates[1:-1]:
            # coupons paid on a settlement date are paid to the seller
            if dt > settlement_date:
                flow_str += ("%12s %12.5f \n" % (dt, flow))

        redemption_amount = face + flow
        flow_str += ("%12s %12.5f \n"
                     % (self._coupon_dates[-1], redemption_amount))

        return flow_str

    ###########################################################################

    def print_bond_payments(self,
                            settlement_date: Date,
                            face: (float) = 100.0):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """

        print(self.bond_payments(settlement_date, face))

    ###########################################################################

    def dirty_price_from_survival_curve(self,
                                        settlement_date: Date,
                                        discount_curve: DiscountCurve,
                                        survival_curve: DiscountCurve,
                                        recovery_rate: float):
        """ Calculate discounted present value of flows assuming default model.
        The survival curve treats the coupons as zero recovery payments while
        the recovery fraction of the par amount is paid at default. For the
        defaulting principal we discretize the time steps using the coupon
        payment times. A finer discretization may handle the time value with
        more accuracy. I reduce any error by averaging period start and period
        end payment present values. """

        f = self._frequency
        c = self._coupon

        pv = 0.0
        prevQ = 1.0
        prevDf = 1.0

        defaultingPrincipalPVPayStart = 0.0
        defaultingPrincipalPVPayEnd = 0.0

        for dt in self._coupon_dates[1:]:

            # coupons paid on a settlement date are paid to the seller
            if dt > settlement_date:
                df = discount_curve.df(dt)
                q = survival_curve.survival_prob(dt)

                # Add PV of coupon conditional on surviving to payment date
                # Any default results in all subsequent coupons being lost
                # with zero recovery

                pv = pv + (c / f) * df * q
                dq = q - prevQ

                defaultingPrincipalPVPayStart += -dq * recovery_rate * prevDf
                defaultingPrincipalPVPayStart += -dq * recovery_rate * df

                # Add on PV of principal if default occurs in coupon period
                prevQ = q
                prevDf = df

        pv = pv + 0.50 * defaultingPrincipalPVPayStart
        pv = pv + 0.50 * defaultingPrincipalPVPayEnd
        pv = pv + df * q * self._redemption
        pv *= self._par
        return pv

    ###########################################################################

    def clean_price_from_survival_curve(self,
                                        settlement_date: Date,
                                        discount_curve: DiscountCurve,
                                        survival_curve: DiscountCurve,
                                        recovery_rate: float):
        """ Calculate clean price value of flows assuming default model.
        The survival curve treats the coupons as zero recovery payments while
        the recovery fraction of the par amount is paid at default. """

        self.accrued_interest(settlement_date, 1.0)

        dirty_price = self.dirty_price_from_survival_curve(settlement_date,
                                                           discount_curve,
                                                           survival_curve,
                                                           recovery_rate)

        clean_price = dirty_price - self._accrued_interest
        return clean_price

    ###########################################################################

    def calc_ror(self,
                 begin_date: Date,
                 end_date: Date,
                 begin_ytm: float,
                 end_ytm: float,
                 convention: YTMCalcType = YTMCalcType.US_STREET):
        """
        Calculate the rate of total return(capital return and interest) given a BUY YTM and a SELL YTM of this bond.
        This function computes the full prices at buying and selling, plus the coupon payments during the period.
        It returns a tuple which includes a simple rate of return, a compounded IRR and the PnL.
        """
        buy_price = self.dirty_price_from_ytm(begin_date, begin_ytm, convention)
        sell_price = self.dirty_price_from_ytm(end_date, end_ytm, convention)
        dates_cfs = zip(self._coupon_dates, self._flow_amounts)
        # The coupon or par payments on buying date belong to the buyer.
        # The coupon or par payments on selling date are given to the new buyer.
        dates_cfs = [(d, c * self._par) for (d, c) in dates_cfs if (d >= begin_date) and (d < end_date)]
        dates_cfs.append((begin_date, -buy_price))
        dates_cfs.append((end_date, sell_price))
        times_cfs = [((d - begin_date)/365, c) for (d, c) in dates_cfs]
        pnl = sum(c for (t, c) in times_cfs)
        simple_return = (pnl / buy_price) * 365 / (end_date - begin_date)
        brentq_up_bound = 5
        brentq_down_bound = -0.9999

        # in case brentq cannot find the irr root
        if simple_return > brentq_up_bound or simple_return < brentq_down_bound:
            irr = simple_return
        else:
            irr = optimize.brentq(npv,
                                  a=brentq_down_bound,  # f(a) and f(b) must have opposite signs
                                  b=brentq_up_bound,
                                  xtol=1e-8,
                                  args=(np.array(times_cfs),)
                                  )

        return simple_return, irr, pnl

    ###########################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("ISSUE DATE", self._issue_date)
        s += label_to_string("MATURITY DATE", self._maturity_date)
        s += label_to_string("COUPON (%)", self._coupon*100.0)
        s += label_to_string("FREQUENCY", self._freq_type)
        s += label_to_string("ACCRUAL TYPE", self._accrual_type)
        s += label_to_string("EX_DIV DAYS", self._ex_div_days, "")
        return s

    ###########################################################################

    def _print(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        print(self)

###############################################################################
