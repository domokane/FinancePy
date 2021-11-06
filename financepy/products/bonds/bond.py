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
from ...market.curves.discount_curve import DiscountCurve


# References https://www.dmo.gov.uk/media/15011/yldeqns_v1.pdf
# DO TRUE YIELD
# JAPANESE SIMPLE YIELD

###############################################################################


from enum import Enum


class YTMCalcType(Enum):
    UK_DMO = 1,
    US_STREET = 2,
    US_TREASURY = 3


###############################################################################


def _f(y, *args):
    """ Function used to do root search in price to yield calculation. """
    bond = args[0]
    settlement_date = args[1]
    price = args[2]
    convention = args[3]
    px = bond.full_price_from_ytm(settlement_date, y, convention)
    obj_fn = px - price
    return obj_fn


###############################################################################


def _g(oas, *args):
    """ Function used to do root search in price to OAS calculation. """
    bond = args[0]
    settlement_date = args[1]
    price = args[2]
    discount_curve = args[3]
    px = bond.full_price_from_oas(settlement_date, discount_curve, oas)
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
                 face_amount: float = 100.0):
        """ Create Bond object by providing the issue date, maturity Date,
        coupon frequency, annualised coupon, the accrual convention type, face
        amount and the number of ex-dividend days. """

        check_argument_types(self.__init__, locals())

        if issue_date >= maturity_date:
            raise FinError("Issue Date must preceded maturity date.")

        self._issue_date = issue_date
        self._maturity_date = maturity_date
        self._coupon = coupon
        self._freq_type = freq_type
        self._accrual_type = accrual_type
        self._frequency = annual_frequency(freq_type)
        self._face_amount = face_amount  # This is the bond holding size
        self._par = 100.0  # This is how price is quoted and amount at maturity
        self._redemption = 1.0  # This is amount paid at maturity

        self._flow_dates = []
        self._flow_amounts = []

        self._accrued_interest = None
        self._accrued_days = 0.0
        self._alpha = 0.0

        self._calculate_flow_dates()
        self._calculate_flows()

    ###########################################################################

    def _calculate_flow_dates(self):
        """ Determine the bond cash flow payment dates."""

        # This should only be called once from init

        calendar_type = CalendarTypes.NONE
        bus_day_rule_type = BusDayAdjustTypes.NONE
        date_gen_rule_type = DateGenRuleTypes.BACKWARD

        self._flow_dates = Schedule(self._issue_date,
                                    self._maturity_date,
                                    self._freq_type,
                                    calendar_type,
                                    bus_day_rule_type,
                                    date_gen_rule_type)._generate()

    ###########################################################################

    def _calculate_flows(self):
        """ Determine the bond cash flow payment amounts without principal """

        self._flow_amounts = [0.0]

        for _ in self._flow_dates[1:]:
            cpn = self._coupon / self._frequency
            self._flow_amounts.append(cpn)

    ###########################################################################

    def full_price_from_ytm(self,
                            settlement_date: Date,
                            ytm: float,
                            convention: YTMCalcType = YTMCalcType.UK_DMO):
        """ Calculate the full price of bond from its yield to maturity. This
        function is vectorised with respect to the yield input. It implements
        a number of standard conventions for calculating the YTM. """

        if convention not in YTMCalcType:
            raise FinError("Yield convention unknown." + str(convention))

        self.calc_accrued_interest(settlement_date)

        ytm = np.array(ytm)  # VECTORIZED
        ytm = ytm + 0.000000000012345  # SNEAKY LOW-COST TRICK TO AVOID y=0

        f = annual_frequency(self._freq_type)
        c = self._coupon
        v = 1.0 / (1.0 + ytm / f)

        # n is the number of flows after the next coupon
        n = 0
        for dt in self._flow_dates:
            if dt > settlement_date:
                n += 1
        n = n - 1

        if n < 0:
            raise FinError("No coupons left")

        if convention == YTMCalcType.UK_DMO:
            if n == 0:
                fp = (v ** (self._alpha)) * (self._redemption + c / f)
            else:
                term1 = (c / f)
                term2 = (c / f) * v
                term3 = (c / f) * v * v * (1.0 - v ** (n - 1)) / (1.0 - v)
                term4 = self._redemption * (v ** n)
                fp = (v ** (self._alpha)) * (term1 + term2 + term3 + term4)
        elif convention == YTMCalcType.US_TREASURY:
            if n == 0:
                fp = (v ** (self._alpha)) * (self._redemption + c / f)
            else:
                term1 = (c / f)
                term2 = (c / f) * v
                term3 = (c / f) * v * v * (1.0 - v ** (n - 1)) / (1.0 - v)
                term4 = self._redemption * (v ** n)
                vw = 1.0 / (1.0 + self._alpha * ytm / f)
                fp = (vw) * (term1 + term2 + term3 + term4)
        elif convention == YTMCalcType.US_STREET:
            vw = 1.0 / (1.0 + self._alpha * ytm / f)
            if n == 0:
                vw = 1.0 / (1.0 + self._alpha * ytm / f)
                fp = vw * (self._redemption + c / f)
            else:
                term1 = (c / f)
                term2 = (c / f) * v
                term3 = (c / f) * v * v * (1.0 - v ** (n - 1)) / (1.0 - v)
                term4 = self._redemption * (v ** n)
                fp = (v ** (self._alpha)) * (term1 + term2 + term3 + term4)
        else:
            raise FinError("Unknown yield convention")

        return fp * self._par

    ###########################################################################

    def principal(self,
                  settlement_date: Date,
                  y: float,
                  convention: YTMCalcType):
        """ Calculate the principal value of the bond based on the face
        amount from its discount margin and making assumptions about the
        future Ibor rates. """

        full_price = self.full_price_from_ytm(settlement_date, y, convention)

        principal = full_price * self._face_amount / self._par
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
        p0 = self.full_price_from_ytm(settlement_date, ytm - dy, convention)
        p2 = self.full_price_from_ytm(settlement_date, ytm + dy, convention)
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
        fp = self.full_price_from_ytm(settlement_date, ytm, convention)
        md = dd * (1.0 + ytm / self._frequency) / fp
        return md

    ###########################################################################

    def modified_duration(self,
                          settlement_date: Date,
                          ytm: float,
                          convention: YTMCalcType = YTMCalcType.UK_DMO):
        """ Calculate the modified duration of the bondon a settlement date
        given its yield to maturity. """

        dd = self.dollar_duration(settlement_date, ytm, convention)
        fp = self.full_price_from_ytm(settlement_date, ytm, convention)
        md = dd / fp
        return md

    ###########################################################################

    def convexity_from_ytm(self,
                           settlement_date: Date,
                           ytm: float,
                           convention: YTMCalcType = YTMCalcType.UK_DMO):
        """ Calculate the bond convexity from the yield to maturity. This
        function is vectorised with respect to the yield input. """

        dy = 0.0001
        p0 = self.full_price_from_ytm(settlement_date, ytm - dy, convention)
        p1 = self.full_price_from_ytm(settlement_date, ytm, convention)
        p2 = self.full_price_from_ytm(settlement_date, ytm + dy, convention)
        conv = ((p2 + p0) - 2.0 * p1) / dy / dy / p1 / self._par
        return conv

    ###########################################################################

    def clean_price_from_ytm(self,
                             settlement_date: Date,
                             ytm: float,
                             convention: YTMCalcType = YTMCalcType.UK_DMO):
        """ Calculate the bond clean price from the yield to maturity. This
        function is vectorised with respect to the yield input. """

        full_price = self.full_price_from_ytm(settlement_date, ytm, convention)
        accrued_amount = self._accrued_interest * self._par / self._face_amount
        clean_price = full_price - accrued_amount
        return clean_price

    ###########################################################################

    def clean_price_from_discount_curve(self,
                                        settlement_date: Date,
                                        discount_curve: DiscountCurve):
        """ Calculate the clean bond value using some discount curve to
        present-value the bond's cash flows back to the curve anchor date and
        not to the settlement date. """

        self.calc_accrued_interest(settlement_date)
        full_price = self.full_price_from_discount_curve(settlement_date,
                                                         discount_curve)

        accrued = self._accrued_interest * self._par / self._face_amount
        clean_price = full_price - accrued
        return clean_price

    ###########################################################################

    def full_price_from_discount_curve(self,
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

        px = 0.0
        df = 1.0
        dfSettle = discount_curve.df(settlement_date)

        for dt in self._flow_dates[1:]:

            # coupons paid on the settlement date are included
            if dt >= settlement_date:
                df = discount_curve.df(dt)
                flow = self._coupon / self._frequency
                pv = flow * df
                px += pv

        px += df * self._redemption
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

        self.calc_accrued_interest(settlement_date)
        accrued_amount = self._accrued_interest * self._par / self._face_amount
        full_prices = (clean_prices + accrued_amount)
        ytms = []

        for full_price in full_prices:
            argtuple = (self, settlement_date, full_price, convention)

            ytm = optimize.newton(_f,
                                  x0=0.05,  # guess initial value of 10%
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

    def calc_accrued_interest(self,
                              settlement_date: Date,
                              num_ex_dividend_days: int = 0,
                              calendar_type: CalendarTypes = CalendarTypes.WEEKEND):
        """ Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date. Note that for some day
        count schemes (such as 30E/360) this is not actually the number of days
        between the previous coupon payment date and settlement date. If the
        bond trades with ex-coupon dates then you need to supply the number of
        days before the coupon date the ex-coupon date is. You can specify the
        calendar to be used - NONE means only calendar days, WEEKEND is only
        weekends or you can specify a country calendar for business days."""

        num_flows = len(self._flow_dates)

        if num_flows == 0:
            raise FinError("Accrued interest - not enough flow dates.")

        for iFlow in range(1, num_flows):
            # coupons paid on a settlement date are paid
            if self._flow_dates[iFlow] >= settlement_date:
                self._pcd = self._flow_dates[iFlow - 1]
                self._ncd = self._flow_dates[iFlow]
                break

        dc = DayCount(self._accrual_type)
        cal = Calendar(calendar_type)
        exDividend_date = cal.add_business_days(
            self._ncd, -num_ex_dividend_days)

        (acc_factor, num, _) = dc.year_frac(self._pcd,
                                            settlement_date,
                                            self._ncd,
                                            self._freq_type)

        if settlement_date > exDividend_date:
            acc_factor = acc_factor - 1.0 / self._frequency

        self._alpha = 1.0 - acc_factor * self._frequency
        self._accrued_interest = acc_factor * self._face_amount * self._coupon
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
        self.calc_accrued_interest(settlement_date)
        accrued_amount = self._accrued_interest * self._par / self._face_amount
        bondPrice = clean_price + accrued_amount
        # Calculate the price of the bond discounted on the Ibor curve
        pvIbor = 0.0
        prev_date = self._pcd

        for dt in self._flow_dates[1:]:

            # coupons paid on a settlement date are included
            if dt >= settlement_date:
                df = discount_curve.df(dt)
                pvIbor += df * self._coupon / self._frequency

        pvIbor += df * self._redemption

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

    def full_price_from_oas(self,
                            settlement_date: Date,
                            discount_curve: DiscountCurve,
                            oas: float):
        """ Calculate the full price of the bond from its OAS given the bond
        settlement date, a discount curve and the oas as a number. """

        self.calc_accrued_interest(settlement_date)
        f = self._frequency
        c = self._coupon

        pv = 0.0
        for dt in self._flow_dates[1:]:

            # coupons paid on a settlement date are included
            if dt >= settlement_date:
                t = (dt - settlement_date) / gDaysInYear

                t = np.maximum(t, gSmall)

                df = discount_curve.df(dt)
                # determine the Ibor implied zero rate
                r = f * (np.power(df, -1.0 / t / f) - 1.0)
                # determine the OAS adjusted zero rate
                df_adjusted = np.power(1.0 + (r + oas) / f, -t * f)
                pv = pv + (c / f) * df_adjusted

        pv = pv + df_adjusted * self._redemption
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

        self.calc_accrued_interest(settlement_date)

        accrued_amount = self._accrued_interest * self._par / self._face_amount
        full_prices = clean_prices + accrued_amount

        oass = []

        for full_price in full_prices:
            argtuple = (self, settlement_date, full_price, discount_curve)

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

    def print_flows(self,
                    settlement_date: Date):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """

        flow = self._face_amount * self._coupon / self._frequency

        for dt in self._flow_dates[1:-1]:
            # coupons paid on a settlement date are included
            if dt >= settlement_date:
                print("%12s" % dt, " %12.2f " % flow)

        redemption_amount = self._face_amount + flow
        print("%12s" % self._flow_dates[-1], " %12.2f " % redemption_amount)

    ###########################################################################

    def full_price_from_survival_curve(self,
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

        for dt in self._flow_dates[1:]:

            # coupons paid on a settlement date are included
            if dt >= settlement_date:
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

        self.calc_accrued_interest(settlement_date)

        full_price = self.full_price_from_survival_curve(settlement_date,
                                                         discount_curve,
                                                         survival_curve,
                                                         recovery_rate)

        clean_price = full_price - self._accrued_interest
        return clean_price

    ###########################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("ISSUE DATE", self._issue_date)
        s += label_to_string("MATURITY DATE", self._maturity_date)
        s += label_to_string("COUPON", self._coupon)
        s += label_to_string("FREQUENCY", self._freq_type)
        s += label_to_string("ACCRUAL TYPE", self._accrual_type)
        s += label_to_string("FACE AMOUNT", self._face_amount, "")
        return s

    ###########################################################################

    def _print(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        print(self)

###############################################################################
