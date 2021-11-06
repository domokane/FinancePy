##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from scipy import optimize

from ...utils.date import Date
from ...utils.error import FinError
from ...utils.frequency import annual_frequency, FrequencyTypes
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.schedule import Schedule
from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes
from ...utils.calendar import DateGenRuleTypes
from ...utils.helpers import label_to_string, check_argument_types


###############################################################################
# TODO: Need to complete and verify the risk sensitivity calculations.
###############################################################################


def _f(dm, *args):
    """ Function used to do solve root search in DM calculation """

    self = args[0]
    settlement_date = args[1]
    next_coupon = args[2]
    current_ibor = args[3]
    future_ibor = args[4]
    full_price = args[5]

    px = self.full_price_from_dm(settlement_date,
                                 next_coupon,
                                 current_ibor,
                                 future_ibor,
                                 dm)

    obj_fn = px - full_price
    return obj_fn


###############################################################################


class BondFRN:
    """ Class for managing floating rate notes that pay a floating index plus a
    quoted margin."""

    def __init__(self,
                 issue_date: Date,
                 maturity_date: Date,
                 quoted_margin: float,  # Fixed spread paid on top of index
                 freq_type: FrequencyTypes,
                 accrual_type: DayCountTypes,
                 face_amount: float = 100.0):
        """ Create FinFloatingRateNote object given its maturity date, its
        quoted margin, coupon frequency, accrual type. Face is the size of
        the position and par is the notional on which price is quoted. """

        check_argument_types(self.__init__, locals())

        self._issue_date = issue_date
        self._maturity_date = maturity_date
        self._quoted_margin = quoted_margin
        self._freq_type = freq_type
        self._accrual_type = accrual_type
        self._flow_dates = []
        self._frequency = annual_frequency(freq_type)
        self._face_amount = face_amount  # This is the position size
        self._par = 100.0  # This is how price is quoted
        self._redemption = 1.0  # This is amount paid at maturity TODO NOT USED

        self._flow_dates = []
        self._flow_amounts = []

        self._settlement_date = Date(1, 1, 1900)
        self._accrued_interest = None
        self._accrued_days = 0.0

        self._calculate_flow_dates()

    ###############################################################################

    def _calculate_flow_dates(self):
        """ Determine the bond cashflow payment dates. """

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

    ###############################################################################

    def full_price_from_dm(self,
                           settlement_date: Date,
                           next_coupon: float,  # The total reset coupon on NCD
                           current_ibor: float,  # Ibor discount to NCD
                           future_ibor: float,  # Future constant Ibor rates
                           dm: float):  # Discount margin
        """ Calculate the full price of the bond from its discount margin (DM)
        using standard model based on assumptions about future Ibor rates. The
        next Ibor payment which has reset is entered, so to is the current
        Ibor rate from settlement to the next coupon date (NCD). Finally there
        is the level of subsequent future Ibor payments and the discount
        margin. """

        self.calc_accrued_interest(settlement_date, next_coupon)

        day_counter = DayCount(self._accrual_type)

        q = self._quoted_margin
        num_flows = len(self._flow_dates)

        # We discount using Libor over the period from settlement to the ncd
        (alpha, _, _) = day_counter.year_frac(settlement_date, self._ncd)
        df = 1.0 / (1.0 + alpha * (current_ibor + dm))

        # A full coupon is paid
        (alpha, _, _) = day_counter.year_frac(self._pcd, self._ncd)
        pv = next_coupon * alpha * df

        # Now do all subsequent coupons that fall after the ncd
        for iFlow in range(1, num_flows):

            if self._flow_dates[iFlow] > self._ncd:
                pcd = self._flow_dates[iFlow - 1]
                ncd = self._flow_dates[iFlow]
                (alpha, _, _) = day_counter.year_frac(pcd, ncd)

                df = df / (1.0 + alpha * (future_ibor + dm))
                c = future_ibor + q
                pv = pv + c * alpha * df

        pv += df
        pv = pv * self._par
        return pv

    ###############################################################################

    def principal(self,
                  settlement_date: Date,
                  next_coupon: float,
                  current_ibor: float,
                  future_ibor: float,
                  dm: float):
        """ Calculate the clean trade price of the bond based on the face
        amount from its discount margin and making assumptions about the
        future Ibor rates. """

        full_price = self.full_price_from_dm(settlement_date,
                                             next_coupon,
                                             current_ibor,
                                             future_ibor,
                                             dm)

        accrued = self._accrued_interest
        principal = full_price * self._face_amount / self._par - accrued
        return principal

    ###############################################################################

    def dollar_duration(self,
                        settlement_date: Date,
                        next_coupon: float,
                        current_ibor: float,
                        future_ibor: float,
                        dm: float):
        """ Calculate the risk or dP/dy of the bond by bumping. This is also
        known as the DV01 in Bloomberg. """

        dy = 0.0001  # 1 basis point

        p0 = self.full_price_from_dm(settlement_date,
                                     next_coupon,
                                     current_ibor + dy,
                                     future_ibor,
                                     dm)

        p2 = self.full_price_from_dm(settlement_date,
                                     next_coupon,
                                     current_ibor - dy,
                                     future_ibor,
                                     dm)

        durn = (p2 - p0) / dy / 2.0
        return durn

    ###############################################################################

    def dollar_credit_duration(self,
                               settlement_date: Date,
                               next_coupon: float,
                               current_ibor: float,
                               future_ibor: float,
                               dm: float):
        """ Calculate the risk or dP/dy of the bond by bumping. """

        if dm > 10.0:
            raise FinError("Discount margin exceeds 100000bp")

        self.calc_accrued_interest(settlement_date, next_coupon)
        dy = 0.0001

        p0 = self.full_price_from_dm(settlement_date,
                                     next_coupon,
                                     current_ibor,
                                     future_ibor,
                                     dm + dy)

        p2 = self.full_price_from_dm(settlement_date,
                                     next_coupon,
                                     current_ibor,
                                     future_ibor,
                                     dm)

        durn = (p2 - p0) / dy
        return durn

    ###############################################################################

    def macauley_duration(self,
                          settlement_date: Date,
                          next_coupon: float,
                          current_ibor: float,
                          future_ibor: float,
                          dm: float):
        """ Calculate the Macauley duration of the FRN on a settlement date
        given its yield to maturity. """

        dd = self.dollar_duration(settlement_date,
                                  next_coupon,
                                  current_ibor,
                                  future_ibor,
                                  dm)

        fp = self.full_price_from_dm(settlement_date,
                                     next_coupon,
                                     current_ibor,
                                     future_ibor,
                                     dm)

        md = dd * (1.0 + (next_coupon + dm) / self._frequency) / fp
        return md

    ###############################################################################

    def modified_duration(self,
                          settlement_date: Date,
                          next_coupon: float,
                          current_ibor: float,
                          future_ibor: float,
                          dm: float):
        """ Calculate the modified duration of the bond on a settlement date
        using standard model based on assumptions about future Ibor rates. The
        next Ibor payment which has reset is entered, so to is the current
        Ibor rate from settlement to the next coupon date (NCD). Finally there
        is the level of subsequent future Ibor payments and the discount
        margin. """

        dd = self.dollar_duration(settlement_date,
                                  next_coupon,
                                  current_ibor,
                                  future_ibor,
                                  dm)

        fp = self.full_price_from_dm(settlement_date,
                                     next_coupon,
                                     current_ibor,
                                     future_ibor,
                                     dm)
        md = dd / fp
        return md

    ###############################################################################

    def modified_credit_duration(self,
                                 settlement_date: Date,
                                 next_coupon: float,
                                 current_ibor: float,
                                 future_ibor: float,
                                 dm: float):
        """ Calculate the modified duration of the bond on a settlement date
        using standard model based on assumptions about future Ibor rates. The
        next Ibor payment which has reset is entered, so to is the current
        Ibor rate from settlement to the next coupon date (NCD). Finally there
        is the level of subsequent future Ibor payments and the discount
        margin. """

        dd = self.dollar_credit_duration(settlement_date,
                                         next_coupon,
                                         current_ibor,
                                         future_ibor,
                                         dm)

        fp = self.full_price_from_dm(settlement_date,
                                     next_coupon,
                                     current_ibor,
                                     future_ibor,
                                     dm)
        md = dd / fp
        return md

    ###############################################################################

    def convexity_from_dm(self,
                          settlement_date: Date,
                          next_coupon: float,
                          current_ibor: float,
                          future_ibor: float,
                          dm: float):
        """ Calculate the bond convexity from the discount margin (DM) using a
        standard model based on assumptions about future Ibor rates. The
        next Ibor payment which has reset is entered, so to is the current
        Ibor rate from settlement to the next coupon date (NCD). Finally there
        is the level of subsequent future Ibor payments and the discount
        margin. """

        dy = 0.0001

        p0 = self.full_price_from_dm(settlement_date,
                                     next_coupon,
                                     current_ibor - dy,
                                     future_ibor,
                                     dm)

        p1 = self.full_price_from_dm(settlement_date,
                                     next_coupon,
                                     current_ibor,
                                     future_ibor,
                                     dm)

        p2 = self.full_price_from_dm(settlement_date,
                                     next_coupon,
                                     current_ibor + dy,
                                     future_ibor,
                                     dm)

        conv = ((p2 + p0) - 2.0 * p1) / dy / dy / p1 / self._par
        return conv

    ###############################################################################

    def clean_price_from_dm(self,
                            settlement_date: Date,
                            next_coupon: float,
                            current_ibor: float,
                            future_ibor: float,
                            dm: float):
        """ Calculate the bond clean price from the discount margin
        using standard model based on assumptions about future Ibor rates. The
        next Ibor payment which has reset is entered, so to is the current
        Ibor rate from settlement to the next coupon date (NCD). Finally there
        is the level of subsequent future Ibor payments and the discount
        margin. """

        if dm > 10.0:
            raise FinError("Discount margin exceeds 100000bp")

        full_price = self.full_price_from_dm(settlement_date,
                                             next_coupon,
                                             current_ibor,
                                             future_ibor,
                                             dm)

        accrued = self._accrued_interest(settlement_date, next_coupon)
        accrued = accrued * self._par / self._face_amount

        clean_price = full_price - accrued
        return clean_price

    ###############################################################################

    def discount_margin(self,
                        settlement_date: Date,
                        next_coupon: float,
                        current_ibor: float,
                        future_ibor: float,
                        clean_price: float):
        """ Calculate the bond's yield to maturity by solving the price
        yield relationship using a one-dimensional root solver. """

        self.calc_accrued_interest(settlement_date, next_coupon)

        # Needs to be adjusted to par notional
        accrued = self._accrued_interest * self._par / self._face_amount

        full_price = clean_price + accrued

        argtuple = (self, settlement_date, next_coupon, current_ibor,
                    future_ibor, full_price)

        dm = optimize.newton(_f,
                             x0=0.01,  # initial value of 10%
                             fprime=None,
                             args=argtuple,
                             tol=1e-12,
                             maxiter=50,
                             fprime2=None)

        return dm

    ###############################################################################

    def calc_accrued_interest(self,
                              settlement_date: Date,
                              next_coupon: float):
        """ Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date. Ex-dividend dates are 
        not handled. Contact me if you need this functionality. """

        num_flows = len(self._flow_dates)

        if num_flows == 0:
            raise FinError("Accrued interest - not enough flow dates.")

        dc = DayCount(self._accrual_type)

        for i in range(1, num_flows):
            if self._flow_dates[i] > settlement_date:
                self._pcd = self._flow_dates[i - 1]
                self._ncd = self._flow_dates[i]
                break

        (acc_factor, num, _) = dc.year_frac(self._pcd,
                                            settlement_date,
                                            self._ncd,
                                            self._freq_type)

        self._alpha = 1.0 - acc_factor * self._frequency

        self._accrued_interest = acc_factor * self._face_amount * next_coupon
        self._accrued_days = num
        return self._accrued_interest

    ###############################################################################

    def print_flows(self,
                    settlement_date: Date):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        self._calculate_flow_dates()
        for dt in self._flow_dates[1:-1]:
            print(dt)

        print(self._flow_dates[-1])

    ###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("ISSUE DATE", self._issue_date)
        s += label_to_string("MATURITY DATE", self._maturity_date)
        s += label_to_string("QUOTED MARGIN (bp)",
                             self._quoted_margin * 10000.0)
        s += label_to_string("FREQUENCY", self._freq_type)
        s += label_to_string("ACCRUAL TYPE", self._accrual_type)
        s += label_to_string("FACE AMOUNT", self._face_amount)
        return s

    ###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)

###############################################################################
