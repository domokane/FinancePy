##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from ...utils.Date import Date
from ...utils.Frequency import Frequency, FinFrequencyTypes
from ...utils.Calendar import FinCalendarTypes
from ...utils.Schedule import Schedule
from ...utils.Calendar import FinBusDayAdjustTypes
from ...utils.Calendar import FinDateGenRuleTypes
from ...utils.DayCount import DayCount, FinDayCountTypes
from ...utils.FinError import FinError
from ...utils.FinHelperFunctions import checkArgumentTypes, labelToString
from ...market.curves.FinDiscountCurve import FinDiscountCurve

###############################################################################


class bond_annuity(object):
    """ An annuity is a vector of dates and flows generated according to ISDA
    standard rules which starts on the next date after the start date
    (effective date) and runs up to an end date with no principal repayment.
    Dates are then adjusted according to a specified calendar. """

    def __init__(self,
                 maturity_date: Date,
                 coupon: float,
                 freq_type: FinFrequencyTypes,
                 calendar_type: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 bus_day_adjust_type: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD,
                 dayCountConventionType: FinDayCountTypes = FinDayCountTypes.ACT_360,
                 face: float = 100.0):

        checkArgumentTypes(self.__init__, locals())

        self._maturity_date = maturity_date
        self._coupon = coupon
        self._freq_type = freq_type
        self._frequency = Frequency(freq_type)

        # ISDA Style conventions
        self._calendar_type = calendar_type
        self._bus_day_adjust_type = bus_day_adjust_type
        self._date_gen_rule_type = date_gen_rule_type
        self._dayCountConventionType = dayCountConventionType

        self._face = face
        self._par = 100.0

        self._flow_dates = []
        self._settlement_date = Date(1, 1, 1900)
        self._accruedInterest = None
        self._accrued_days = 0.0
        self._alpha = 0.0

###############################################################################

    def cleanPriceFromDiscountCurve(self,
                                    settlement_date: Date,
                                    discount_curve: FinDiscountCurve):
        """ Calculate the bond price using some discount curve to present-value
        the bond's cashflows. """

        fullPrice = self.fullPriceFromDiscountCurve(settlement_date,
                                                    discount_curve)
        accrued = self._accruedInterest * self._par / self._face
        cleanPrice = fullPrice - accrued
        return cleanPrice

###############################################################################

    def fullPriceFromDiscountCurve(self,
                                   settlement_date: Date,
                                   discount_curve: FinDiscountCurve):
        """ Calculate the bond price using some discount curve to present-value
        the bond's cashflows. """

        self.calculateFlowDatesPayments(settlement_date)
        pv = 0.0

        numFlows = len(self._flow_dates)

        for i in range(1, numFlows):
            dt = self._flow_dates[i]
            df = discount_curve.df(dt)
            flow = self._flow_amounts[i]
            pv = pv + flow * df

        return pv * self._par / self._face

###############################################################################

    def calculateFlowDatesPayments(self,
                                   settlement_date: Date):

        # No need to generate flows if settlement date has not changed
        if settlement_date == self._settlement_date:
            return

        if settlement_date == self._maturity_date:
            raise FinError("Settlement date is maturity date.")

        self._settlement_date = settlement_date
        calendar_type = FinCalendarTypes.NONE
        busDayRuleType = FinBusDayAdjustTypes.NONE
        date_gen_rule_type = FinDateGenRuleTypes.BACKWARD

        self._flow_dates = Schedule(settlement_date,
                                   self._maturity_date,
                                   self._freq_type,
                                   calendar_type,
                                   busDayRuleType,
                                   date_gen_rule_type)._generate()

        self._pcd = self._flow_dates[0]
        self._ncd = self._flow_dates[1]
        self.calcAccruedInterest(settlement_date)

        self._flow_amounts = [0.0]
        basis = DayCount(self._dayCountConventionType)

        prevDt = self._pcd

        for nextDt in self._flow_dates[1:]:
            alpha = basis.year_frac(prevDt, nextDt)[0]
            flow = self._coupon * alpha * self._face
            self._flow_amounts.append(flow)
            prevDt = nextDt

###############################################################################

    def calcAccruedInterest(self,
                            settlement_date: Date):
        """ Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date. """

        if settlement_date != self._settlement_date:
            self.calculateFlowDatesPayments(settlement_date)

        if len(self._flow_dates) == 0:
            raise FinError("Accrued interest - not enough flow dates.")

        dc = DayCount(self._dayCountConventionType)

        (acc_factor, num, _) = dc.year_frac(self._pcd,
                                          settlement_date,
                                          self._ncd, 
                                          self._frequency)

        self._alpha = 1.0 - acc_factor * self._frequency

        self._accruedInterest = acc_factor * self._face * self._coupon
        self._accrued_days = num
        return self._accruedInterest

###############################################################################

    def printFlows(self,
                   settlement_date: Date):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """

        self.calculateFlowDatesPayments(settlement_date)

        numFlows = len(self._flow_dates)
        for i in range(1, numFlows):
            dt = self._flow_dates[i]
            flow = self._flow_amounts[i]
            print(dt, ",", flow)

###############################################################################

    def __repr__(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """

        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("MATURITY DATE", self._maturity_date)
        s += labelToString("FREQUENCY", self._freq_type)
        s += labelToString("CALENDAR", self._calendar_type)
        s += labelToString("BUS_DAY_RULE", self._bus_day_adjust_type)
        s += labelToString("DATE_GEN_RULE", self._date_gen_rule_type)

        return s

###############################################################################

    def _print(self):
        """ Simple print function for backward compatibility. """
        print(self)


###############################################################################
