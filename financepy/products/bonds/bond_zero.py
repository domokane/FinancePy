import numpy as np
from scipy import optimize

from ...utils.date import Date
from ...utils.error import FinError
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
from ...utils.frequency import FrequencyTypes, annual_frequency
from ...products.bonds.bond import YTMCalcType

###############################################################################
# TO DO - THIS CLASS NEEDS TO INHERIT FROM BOND CLASS
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

class BondZero:
    """ A zero coupon bond is a bond which doesn't pay any periodic payments. 
    Instead, it is issued at a discount. The entire face value of the bond is 
    paid out at maturity. It is issued as a deep discount bond. 
    
    There is a special convention for accrued interest in which 
    
        Accrued_interest = (par - issue price) * D
    
    where D = (settlement_date - issue_date)/(maturity_date - issue_date).    
    """

    def __init__(self,
                 issue_date: Date,
                 maturity_date: Date,
                 issue_price: float          # Issue price usually discounted
                 ):
        """ Create BondZero object by providing the issue date, maturity Date,
        face amount and issue price. """

        check_argument_types(self.__init__, locals())

        if issue_date >= maturity_date:
            raise FinError("Issue Date must preceded maturity date.")

        self._issue_date = issue_date
        self._maturity_date = maturity_date
        self._accrual_type = DayCountTypes.ZERO
        self._issue_price = issue_price  # Price of issue, usually discounted
        self._par = 100.0  # This is how price is quoted and amount at maturity
        self._freq_type = FrequencyTypes.ZERO
        self._coupon_dates = [issue_date, maturity_date]
        self._payment_dates = [issue_date, maturity_date]
        self._flow_amounts = [0.0, 0.0]  # coupon payments are zero
        self._calendar_type = CalendarTypes.WEEKEND
        self._ex_div_days = 0
        
        self._accrued_interest = None
        self._accrued_days = 0.0
        self._alpha = 0.0

    ###########################################################################

    def dirty_price_from_ytm(self,
                            settlement_date: Date,
                            ytm: float,
                            convention: YTMCalcType = YTMCalcType.ZERO):
        """ Calculate the full price of bond from its yield to maturity. This
        function is vectorised with respect to the yield input. It implements
        a number of standard conventions for calculating the YTM. """

        if convention != YTMCalcType.ZERO:
            raise FinError("Need to use YTMCalcType.ZERO for zero coupon bond")

        self.accrued_interest(settlement_date, 1.0)

        ytm = np.array(ytm)  # VECTORIZED
        ytm = ytm + 0.000000000012345  # SNEAKY LOW-COST TRICK TO AVOID y=0

        # n is the number of flows after the next coupon
        n = 0
        for dt in self._coupon_dates:
            if dt > settlement_date:
                n += 1
        n = n - 1

        if n < 0:
            raise FinError("No coupons left")
        # A zero coupon bond has a price equal to the discounted principal
        # assuming an annualised rate raised to the power of years

        dc = DayCount(self._accrual_type)
        (acc_factor, num, _) = dc.year_frac(settlement_date,
                                            self._maturity_date,
                                            self._maturity_date,
                                            FrequencyTypes.ZERO)
        if acc_factor <= 1:
            pv = self._par / (1.0 + ytm * acc_factor)
        else:
            pv = self._par / (1.0 + ytm) ** acc_factor
        
        return pv

    ###########################################################################

    def principal(self,
                  settlement_date: Date,
                  ytm: float,
                  face: (float),
                  convention: YTMCalcType = YTMCalcType.ZERO):
        """ Calculate the principal value of the bond based on the face
        amount from its discount margin and making assumptions about the
        future Ibor rates. """

        dirty_price = self.dirty_price_from_ytm(settlement_date, y, convention)

        principal = dirty_price * face / self._par
        principal = principal - self._accrued_interest
        return principal

    ###########################################################################

    def dollar_duration(self,
                        settlement_date: Date,
                        ytm: float,
                        convention: YTMCalcType = YTMCalcType.ZERO):
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
                          convention: YTMCalcType = YTMCalcType.ZERO):
        """ Calculate the Macauley duration of the bond on a settlement date
        given its yield to maturity. """

        dd = self.dollar_duration(settlement_date, ytm, convention)
        fp = self.dirty_price_from_ytm(settlement_date, ytm, convention)
        md = dd * (1.0 + ytm) / fp
        return md

    ###########################################################################

    def modified_duration(self,
                          settlement_date: Date,
                          ytm: float,
                          convention: YTMCalcType = YTMCalcType.ZERO):
        """ Calculate the modified duration of the bondon a settlement date
        given its yield to maturity. """

        dd = self.dollar_duration(settlement_date, ytm, convention)
        fp = self.dirty_price_from_ytm(settlement_date, ytm, convention)
        md = dd / fp
        return md

    ###########################################################################

    def convexity_from_ytm(self,
                           settlement_date: Date,
                           ytm: float,
                           convention: YTMCalcType = YTMCalcType.ZERO):
        """ Calculate the bond convexity from the yield to maturity. This
        function is vectorised with respect to the yield input. """

        dy = 0.0001
        p0 = self.dirty_price_from_ytm(settlement_date, ytm - dy, convention)
        p1 = self.dirty_price_from_ytm(settlement_date, ytm, convention)
        p2 = self.dirty_price_from_ytm(settlement_date, ytm + dy, convention)
        conv = ((p2 + p0) - 2.0 * p1) / dy / dy / p1 / self._par
        return conv

    ###########################################################################

    def clean_price_from_ytm(self,
                             settlement_date: Date,
                             ytm: float,
                             convention: YTMCalcType = YTMCalcType.ZERO):
        """ Calculate the bond clean price from the yield to maturity. This
        function is vectorised with respect to the yield input. """

        dirty_price = self.dirty_price_from_ytm(settlement_date, ytm, convention)
        accrued = self.accrued_interest(settlement_date, self._par)
        clean_price = dirty_price - accrued
        return clean_price

    ###########################################################################

    def clean_price_from_discount_curve(self,
                                        settlement_date: Date,
                                        discount_curve: DiscountCurve):
        """ Calculate the clean bond value using some discount curve to
        present-value the bond's cash flows back to the curve anchor date and
        not to the settlement date. """
 
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

        px = 0.0
        df = 1.0
        dfSettle = discount_curve.df(settlement_date, )

        for dt in self._coupon_dates[1:]:

            # coupons paid on the settlement date are paid to the seller
            if dt > settlement_date:
                df = discount_curve.df(dt)
                flow = 0
                pv = flow * df
                px += pv

        px += df * self._redemption
        px = px / dfSettle

        return px * self._par

    ###########################################################################

    def current_yield(self, clean_price):
        """
        Calculate the current yield of the bond which is the
        coupon divided by the clean price (not the full price). The coupon of a zero coupon bond is defined as:
        (par - issue_price) / tenor
        """
        dc = DayCount(self._accrual_type)
        tenor, _, _ = dc.year_frac(self._issue_date,
                                   self._maturity_date,
                                   self._maturity_date,
                                   FrequencyTypes.ZERO)
        virtual_coupon = (self._par - self._issue_price) / tenor
        y = virtual_coupon / clean_price
        return y

    ###########################################################################

    def yield_to_maturity(self,
                          settlement_date: Date,
                          clean_price: float,
                          convention: YTMCalcType = YTMCalcType.ZERO):
        """ Calculate the bond's yield to maturity by solving the price
        yield relationship using a one-dimensional root solver. """

        if type(clean_price) is float or type(clean_price) is np.float64:
            clean_prices = np.array([clean_price])
        elif type(clean_price) is list or type(clean_price) is np.ndarray:
            clean_prices = np.array(clean_price)
        else:
            raise FinError("Unknown type for clean_price "
                           + str(type(clean_price)))

        accrued_amount = self.accrued_interest(settlement_date, self._par)
        
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

    def accrued_interest(self,
                         settlement_date: Date,
                         face: (float)):
        """ Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date. Note that for some day
        count schemes (such as 30E/360) this is not actually the number of days
        between the previous coupon payment date and settlement date. If the
        bond trades with ex-coupon dates then you need to supply the number of
        days before the coupon date the ex-coupon date is. You can specify the
        calendar to be used - NONE means only calendar days, WEEKEND is only
        weekends or you can specify a country calendar for business days."""

        num_flows = len(self._coupon_dates)

        if num_flows == 0:
            raise FinError("Accrued interest - not enough flow dates.")

        for iFlow in range(1, num_flows):
            # coupons paid on the settlement date are paid to the seller
            if self._coupon_dates[iFlow] > settlement_date:
                self._pcd = self._coupon_dates[iFlow - 1]
                self._ncd = self._coupon_dates[iFlow]
                break

        dc = DayCount(self._accrual_type)
        cal = Calendar(self._calendar_type)
        exDividend_date = cal.add_business_days(
            self._ncd, -self._ex_div_days)

        (acc_factor, num, _) = dc.year_frac(self._pcd,
                                            settlement_date,
                                            self._ncd,
                                            FrequencyTypes.ZERO)

        if settlement_date > exDividend_date:
            acc_factor = acc_factor - 1.0

        self._alpha = 1.0 - acc_factor

        num = (settlement_date - self._issue_date)
        den = (self._maturity_date - self._issue_date)

        f = num / den        
        g = ((self._par - self._issue_price)) / self._par
        
        self._accrued_interest = f * g * face
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

            # coupons paid on the settlement date are paid to the seller
            if dt > settlement_date:
                df = discount_curve.df(dt)
                # pvIbor += df * self._coupon / self._frequency

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

    def dirty_price_from_oas(self,
                            settlement_date: Date,
                            discount_curve: DiscountCurve,
                            oas: float):
        """ Calculate the full price of the bond from its OAS given the bond
        settlement date, a discount curve and the oas as a number. """

        self.accrued_interest(settlement_date, 1.0)

        pv = 0.0
        for dt in self._coupon_dates[1:]:

            # coupons paid on the settlement date are paid to the seller
            if dt > settlement_date:
                t = (dt - settlement_date) / gDaysInYear

                t = np.maximum(t, gSmall)

                df = discount_curve.df(dt)
                # determine the Ibor implied zero rate
                r = np.power(df, -1.0 / t) - 1.0
                # determine the OAS adjusted zero rate
                df_adjusted = np.power(1.0 + (r + oas), -t)

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

    def bond_payments(self,
                      settlement_date: Date, 
                      face: (float)):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        flow_str = ''
        flow_str += ("%12s %12.2f \n"
                     % (self._coupon_dates[-1], face))

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

        pv = 0.0
        prevQ = 1.0
        prevDf = 1.0

        defaultingPrincipalPVPayStart = 0.0
        defaultingPrincipalPVPayEnd = 0.0

        for dt in self._coupon_dates[1:]:

            # coupons paid on the settlement date are paid to the seller
            if dt > settlement_date:
                df = discount_curve.df(dt)
                q = survival_curve.survival_prob(dt)

                # Add PV of coupon conditional on surviving to payment date
                # Any default results in all subsequent coupons being lost
                # with zero recovery

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
                       convention: YTMCalcType = YTMCalcType.ZERO):
        """
        Calculates the rate of total return (capital return and interest) given 
        a BUY YTM and a SELL YTM of this bond.
        
        This function computes the full prices at buying and selling, plus the 
        coupon payments during the period.
        
        It returns a tuple which includes a simple rate of return, a compounded 
        IRR and the PnL.
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

##############################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("ISSUE DATE", self._issue_date)
        s += label_to_string("MATURITY DATE", self._maturity_date)
        s += label_to_string("COUPON (%)", 0)
        s += label_to_string("ISSUE PRICE", self._issue_price)
        s += label_to_string("FREQUENCY", self._freq_type)
        s += label_to_string("ACCRUAL TYPE", self._accrual_type)
        return s

    ###########################################################################

    def _print(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        print(self)

###############################################################################
