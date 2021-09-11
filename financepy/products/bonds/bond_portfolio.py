###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

# UNDER CONSTRUCTION !!!!!!!!!!!!!!!!

import numpy as np

from ...utils.date import Date
from ...utils.calendar import CalendarTypes
from ...utils.helpers import label_to_string, check_argument_types
from ...market.curves.discount_curve import DiscountCurve
from .bond import YTMCalcType

###############################################################################


class BondPortfolio:
    """ Class for valuing and risk-managing a portfolio of bonds. """

    def __init__(self,
                 bonds: (list),
                 bondWeights: (list, np.ndarray)):
        """ XXX """

        check_argument_types(self.__init__, locals())

        self.calculateFlows()

###############################################################################

    def _calculate_flows(self):
        """ Determine the bond cashflow payment amounts without principal """

        self._flow_amounts = [0.0]

        for _ in self._flow_dates[1:]:
            cpn = self._coupon / self._frequency
            self._flow_amounts.append(cpn)

###############################################################################

    def dollar_duration(self,
                        settlement_date: Date,
                        ytm: float,
                        convention: YTMCalcType = YTMCalcType.UK_DMO):
        """ Calculate the risk or dP/dy of the bond by bumping. This is also
        known as the DV01 in Bloomberg. """

        dy = 0.0001  # 1 basis point
        p0 = self.full_priceFromYTM(settlement_date, ytm - dy, convention)
        p2 = self.full_priceFromYTM(settlement_date, ytm + dy, convention)
        durn = -(p2 - p0) / dy / 2.0
        return durn

###############################################################################

    def macauley_duration(self,
                          settlement_date: Date,
                          ytm: float,
                          convention: YTMCalcType = YTMCalcType.UK_DMO):
        """ Calculate the Macauley duration of the bond on a settlement date
        given its yield to maturity. """

        dd = self.dollar_duration(settlement_date, ytm, convention)
        fp = self.full_priceFromYTM(settlement_date, ytm, convention)
        md = dd * (1.0 + ytm / self._frequency) / fp
        return md

###############################################################################

    def modified_duration(self,
                          settlement_date: Date,
                          ytm: float,
                          convention: YTMCalcType = YTMCalcType.UK_DMO):
        """ Calculate the modified duration of the bondon a settlement date
        given its yield to maturity. """

        dd = self.dollar_duration(settlement_date, ytm, convention)
        fp = self.full_priceFromYTM(settlement_date, ytm, convention)
        md = dd / fp
        return md

###############################################################################

    def convexity_from_ytm(self,
                           settlement_date: Date,
                           ytm: float,
                           convention: YTMCalcType = YTMCalcType.UK_DMO):
        """ Calculate the bond convexity from the yield to maturity. This
        function is vectorised with respect to the yield input. """

        dy = 0.0001
        p0 = self.full_priceFromYTM(settlement_date, ytm - dy, convention)
        p1 = self.full_priceFromYTM(settlement_date, ytm, convention)
        p2 = self.full_priceFromYTM(settlement_date, ytm + dy, convention)
        conv = ((p2 + p0) - 2.0 * p1) / dy / dy / p1 / self._par
        return conv

###############################################################################

    def clean_price_from_ytm(self,
                             settlement_date: Date,
                             ytm: float,
                             convention: YTMCalcType = YTMCalcType.UK_DMO):
        """ Calculate the bond clean price from the yield to maturity. This
        function is vectorised with respect to the yield input. """

        full_price = self.full_priceFromYTM(settlement_date, ytm, convention)
        accrued_amount = self._accrued_interest * self._par / self._face_amount
        clean_price = full_price - accrued_amount
        return clean_price

###############################################################################

    def clean_price_from_discount_curve(self,
                                        settlement_date: Date,
                                        discount_curve: DiscountCurve):
        """ Calculate the clean bond value using some discount curve to
        present-value the bond's cash flows back to the curve anchor date and
        not to the settlement date. """


##############################################################################

    def full_price_from_discount_curve(self,
                                       settlement_date: Date,
                                       discount_curve: DiscountCurve):
        """ Calculate the bond price using a provided discount curve to PV the
        bond's cash flows to the settlement date. As such it is effectively a
        forward bond price if the settlement date is after the valuation date.
        """


###############################################################################

    def current_yield(self, clean_price):
        """ Calculate the current yield of the bond which is the
        coupon divided by the clean price (not the full price)"""

        y = self._coupon * self._par / clean_price
        return y

###############################################################################

    def yield_to_maturity(self,
                          settlement_date: Date,
                          clean_price: float,
                          convention: YTMCalcType = YTMCalcType.US_TREASURY):
        """ Calculate the bond's yield to maturity by solving the price
        yield relationship using a one-dimensional root solver. """


###############################################################################

    def calc_accrued_interest(self,
                              settlement_date: Date,
                              num_ex_dividend_days: int = 0,
                              calendar_type: CalendarTypes = CalendarTypes.WEEKEND):

        return self._accrued_interest

###############################################################################

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

###############################################################################

    def full_price_from_survival_curve(self,
                                       settlement_date: Date,
                                       discount_curve: DiscountCurve,
                                       survival_curve: DiscountCurve,
                                       recovery_rate: float):
        """ Calculate discounted present value of flows assuming default model.
        The survival curve treats the coupons as zero recovery payments while
        the recovery fraction of the par amount is paid at default. For the 
        defaulting principal we discretise the time steps using the coupon
        payment times. A finer discretisation may handle the time value with
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

###############################################################################

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

###############################################################################

    def __repr__(self):

        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("ISSUE DATE", self._issue_date)
        s += label_to_string("MATURITY DATE", self._maturity_date)
        s += label_to_string("COUPON", self._coupon)
        s += label_to_string("FREQUENCY", self._freq_type)
        s += label_to_string("ACCRUAL TYPE", self._accrual_type)
        s += label_to_string("FACE AMOUNT", self._face_amount, "")
        return s

###############################################################################

    def _print(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        print(self)


###############################################################################
