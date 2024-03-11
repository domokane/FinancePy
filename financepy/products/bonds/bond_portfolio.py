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


# class BondPortfolio:
#     """ Class for valuing and risk-managing a portfolio of bonds. """

#     def __init__(self,
#                  bonds: (list),
#                  bondWeights: (list, np.ndarray)):
#         """ XXX """

#         check_argument_types(self.__init__, locals())

#         self.calculateFlows()
#         self.par = 100.0

# ###############################################################################

#     def _calculate_flows(self):
#         """ Determine the bond cashflow payment amounts without principal """

#         self.flow_amounts = [0.0]

#         for _ in self.cpn_dts[1:]:
#             cpn = self.cpn / self.freq
#             self.flow_amounts.append(cpn)

# ###############################################################################

#     def dollar_duration(self,
#                         settle_dt: Date,
#                         ytm: float,
#                         convention: YTMCalcType = YTMCalcType.UK_DMO):
#         """ Calculate the risk or dP/dy of the bond by bumping. This is also
#         known as the DV01 in Bloomberg. """

#         dy = 0.0001  # 1 basis point
#         p0 = self.dirty_price_from_ytm(settle_dt, ytm - dy, convention)
#         p2 = self.dirty_price_from_ytm(settle_dt, ytm + dy, convention)
#         durn = -(p2 - p0) / dy / 2.0
#         return durn

# ###############################################################################

#     def macauley_duration(self,
#                           settle_dt: Date,
#                           ytm: float,
#                           convention: YTMCalcType = YTMCalcType.UK_DMO):
#         """ Calculate the Macauley duration of the bond on a settlement date
#         given its yield to maturity. """

#         dd = self.dollar_duration(settle_dt, ytm, convention)
#         fp = self.dirty_price_from_ytm(settle_dt, ytm, convention)
#         md = dd * (1.0 + ytm / self.freq) / fp
#         return md

# ###############################################################################

#     def modified_duration(self,
#                           settle_dt: Date,
#                           ytm: float,
#                           convention: YTMCalcType = YTMCalcType.UK_DMO):
#         """ Calculate the modified duration of the bondon a settlement date
#         given its yield to maturity. """

#         dd = self.dollar_duration(settle_dt, ytm, convention)
#         fp = self.dirty_price_from_ytm(settle_dt, ytm, convention)
#         md = dd / fp
#         return md

# ###############################################################################

#     def convexity_from_ytm(self,
#                            settle_dt: Date,
#                            ytm: float,
#                            convention: YTMCalcType = YTMCalcType.UK_DMO):
#         """ Calculate the bond convexity from the yield to maturity. This
#         function is vectorised with respect to the yield input. """

#         dy = 0.0001
#         p0 = self.dirty_price_from_ytm(settle_dt, ytm - dy, convention)
#         p1 = self.dirty_price_from_ytm(settle_dt, ytm, convention)
#         p2 = self.dirty_price_from_ytm(settle_dt, ytm + dy, convention)
#         conv = ((p2 + p0) - 2.0 * p1) / dy / dy / p1 / self.par
#         return conv

# ###############################################################################

#     def clean_price_from_ytm(self,
#                              settle_dt: Date,
#                              ytm: float,
#                              convention: YTMCalcType = YTMCalcType.UK_DMO):
#         """ Calculate the bond clean price from the yield to maturity. This
#         function is vectorised with respect to the yield input. """

#         dirty_price = self.dirty_price_from_ytm(settle_dt, ytm,
#                                                 convention)
#         accrued_amount = self.accrued_int* self.par
#         clean_price = dirty_price - accrued_amount
#         return clean_price

# ###############################################################################

#     def clean_price_from_discount_curve(self,
#                                         settle_dt: Date,
#                                         discount_curve: DiscountCurve):
#         """ Calculate the clean bond value using some discount curve to
#         present-value the bond's cash flows back to the curve anchor date and
#         not to the settlement date. """


# ##############################################################################

#     def dirty_price_from_discount_curve(self,
#                                         settle_dt: Date,
#                                         discount_curve: DiscountCurve):
#         """ Calculate the bond price using a provided discount curve to PV the
#         bond's cash flows to the settlement date. As such it is effectively a
#         forward bond price if the settlement date is after the valuation date.
#         """


# ###############################################################################

#     def current_yield(self, clean_price):
#         """ Calculate the current yield of the bond which is the
#         coupon divided by the clean price (not the full price)"""

#         y = self.cpn * self.par / clean_price
#         return y

# ###############################################################################

#     def yield_to_maturity(self,
#                           settle_dt: Date,
#                           clean_price: float,
#                           convention: YTMCalcType = YTMCalcType.US_TREASURY):
#         """ Calculate the bond's yield to maturity by solving the price
#         yield relationship using a one-dimensional root solver. """


# ###############################################################################

#     def accrued_interest(self,
#                          settle_dt: Date,
#                          num_ex_dividend_days: int = 0,
#                          cal_type: CalendarTypes = CalendarTypes.WEEKEND):

#         return self.accrued_interest

# ###############################################################################

#     def print_payments(self,
#                        settle_dt: Date,
#                        face: (float)):
#         """ Print a list of the unadjusted coupon payment dates used in
#         analytic calculations for the bond. """

#         flow = self.cpn / self.freq

#         for dt in self.cpn_dts[1:-1]:
#             # coupons paid on a settlement date are included
#             if dt >= settle_dt:
#                 print("%12s" % dt, " %12.2f " % flow)

#         redemption_amount = face * (1.0 + flow)

#         print("%12s" % self.cpn_dts[-1], " %12.2f " % redemption_amount)

# ###############################################################################

#     def dirty_price_from_survival_curve(self,
#                                         settle_dt: Date,
#                                         discount_curve: DiscountCurve,
#                                         survival_curve: DiscountCurve,
#                                         recovery_rate: float):
#         """ Calculate discounted present value of flows assuming default model.
#         The survival curve treats the coupons as zero recovery payments while
#         the recovery fraction of the par amount is paid at default. For the
#         defaulting principal we discretise the time steps using the coupon
#         payment times. A finer discretisation may handle the time value with
#         more accuracy. I reduce any error by averaging period start and period
#         end payment present values. """

#         f = self.freq
#         c = self.cpn

#         pv = 0.0
#         prevQ = 1.0
#         prevDf = 1.0

#         defaultingPrincipalPVPayStart = 0.0
#         defaultingPrincipalPVPayEnd = 0.0

#         for dt in self.cpn_dts[1:]:

#             # coupons paid on a settlement date are included
#             if dt >= settle_dt:

#                 df = discount_curve.df(dt)
#                 q = survival_curve.survival_prob(dt)

#                 # Add PV of coupon conditional on surviving to payment date
#                 # Any default results in all subsequent coupons being lost
#                 # with zero recovery

#                 pv = pv + (c / f) * df * q
#                 dq = q - prevQ

#                 defaultingPrincipalPVPayStart += -dq * recovery_rate * prevDf
#                 defaultingPrincipalPVPayStart += -dq * recovery_rate * df

#                 # Add on PV of principal if default occurs in coupon period
#                 prevQ = q
#                 prevDf = df

#         pv = pv + 0.50 * defaultingPrincipalPVPayStart
#         pv = pv + 0.50 * defaultingPrincipalPVPayEnd
#         pv = pv + df * q * self.par
#         pv *= self.par
#         return pv

# ###############################################################################

#     def clean_price_from_survival_curve(self,
#                                         settle_dt: Date,
#                                         discount_curve: DiscountCurve,
#                                         survival_curve: DiscountCurve,
#                                         recovery_rate: float):
#         """ Calculate clean price value of flows assuming default model.
#         The survival curve treats the coupons as zero recovery payments while
#         the recovery fraction of the par amount is paid at default. """

#         self.accrued_interest(settle_dt, 1.0)

#         dirty_price = self.dirty_price_from_survival_curve(settle_dt,
#                                                            discount_curve,
#                                                            survival_curve,
#                                                            recovery_rate)

#         clean_price = dirty_price - self.accrued_interest
#         return clean_price

# ###############################################################################

#     def __repr__(self):

#         s = label_to_string("OBJECT TYPE", type(self).__name__)
#         s += label_to_string("ISSUE DATE", self.issue_dt)
#         s += label_to_string("MATURITY DATE", self.maturity_dt)
#         s += label_to_string("COUPON", self.cpn)
#         s += label_to_string("FREQUENCY", self.freq_type)
#         s += label_to_string("DAY COUNT TYPE", self.dc_type)
#         return s

# ###############################################################################

#     def _print(self):
#         """ Print a list of the unadjusted coupon payment dates used in
#         analytic calculations for the bond. """
#         print(self)


# ###############################################################################
