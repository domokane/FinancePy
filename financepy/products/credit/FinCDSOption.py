##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from math import sqrt, log
from scipy import optimize


from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes, DateGenRuleTypes
from ...utils.day_count import DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.global_variables import gDaysInYear
from ...utils.fin_math import ONE_MILLION, N
from ...products.credit.cds import FinCDS
from ...utils.helper_functions import check_argument_types
from ...utils.date import Date
from ...utils.FinError import FinError

##########################################################################


def fvol(volatility, *args):
    """ Root searching function in the calculation of the CDS implied
    volatility. """

    self = args[0]
    valuation_date = args[1]
    issuer_curve = args[2]
    optionPrice = args[3]
    value = self.value(valuation_date, issuer_curve, volatility)
    objFn = value - optionPrice

    return objFn

##########################################################################
##########################################################################


class FinCDSOption():
    """ Class to manage the pricing and risk-management of an option on a
    single-name CDS. This is a contract in which the option buyer pays for an
    option to either buy or sell protection on the underlying CDS at a fixed
    spread agreed today and to be exercised in the future on a specified expiry
    date. The option may or may not cancel if there is a credit event before
    option expiry. This needs to be specified. """

    def __init__(self,
                 expiry_date: Date,
                 maturity_date: Date,
                 strikeCoupon: float,
                 notional: float = ONE_MILLION,
                 long_protection: bool = True,
                 knockoutFlag: bool = True,
                 freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 day_count_type: DayCountTypes = DayCountTypes.ACT_360,
                 calendar_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bus_day_adjust_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 date_gen_rule_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):
        """ Create a FinCDSOption object with the option expiry date, the
        maturity date of the underlying CDS, the option strike coupon,
        notional, whether the option knocks out or not in the event of a credit
        event before expiry and the payment details of the underlying CDS. """

        check_argument_types(self.__init__, locals())

        if maturity_date < expiry_date:
            raise FinError("Maturity date must be after option expiry date")

        if strikeCoupon < 0.0:
            raise FinError("Strike must be greater than zero")

        self._expiry_date = expiry_date
        self._maturity_date = maturity_date
        self._strikeCoupon = strikeCoupon
        self._long_protection = long_protection
        self._knockoutFlag = knockoutFlag
        self._notional = notional

        self._freq_type = freq_type
        self._day_count_type = day_count_type
        self._calendar_type = calendar_type
        self._businessDateAdjustType = bus_day_adjust_type
        self._date_gen_rule_type = date_gen_rule_type

##########################################################################

    def value(self,
              valuation_date,
              issuer_curve,
              volatility):
        """ Value the CDS option using Black's model with an adjustment for any
        Front End Protection.
        TODO - Should the CDS be created in the init method ? """

        if valuation_date > self._expiry_date:
            raise FinError("Expiry date is now or in the past")

        if volatility < 0.0:
            raise FinError("Volatility must be greater than zero")

        # The underlying is a forward starting option that steps in on
        # the expiry date and matures on the expiry date with a coupon
        # set equal to the option spread strike
        cds = FinCDS(self._expiry_date,
                     self._maturity_date,
                     self._strikeCoupon,
                     self._notional,
                     self._long_protection,
                     self._freq_type,
                     self._day_count_type,
                     self._calendar_type,
                     self._businessDateAdjustType,
                     self._date_gen_rule_type)

        strike = self._strikeCoupon
        forwardSpread = cds.parSpread(valuation_date, issuer_curve)
        forwardRPV01 = cds.riskyPV01(valuation_date, issuer_curve)['full_rpv01']

        timeToExpiry = (self._expiry_date - valuation_date) / gDaysInYear
        logMoneyness = log(forwardSpread / strike)

        halfVolSquaredT = 0.5 * volatility * volatility * timeToExpiry
        volSqrtT = volatility * sqrt(timeToExpiry)

        d1 = (logMoneyness + halfVolSquaredT) / volSqrtT
        d2 = (logMoneyness - halfVolSquaredT) / volSqrtT

        if self._long_protection:
            optionValue = forwardSpread * N(d1) - strike * N(d2)
        else:
            optionValue = strike * N(-d2) - forwardSpread * N(-d1)

        optionValue = optionValue * forwardRPV01

        # If the option does not knockout on a default before expiry then we
        # need to include the cost of protection which is provided between
        # the value date and the expiry date
        if self._knockoutFlag is False and self._long_protection is True:

            df = issuer_curve.getDF(timeToExpiry)
            q = issuer_curve.getSurvProb(timeToExpiry)
            recovery = issuer_curve._recovery_rate
            frontEndProtection = df * (1.0 - q) * (1.0 - recovery)
            optionValue += frontEndProtection

        # we return the option price in dollars
        return optionValue * self._notional

##########################################################################

    def impliedVolatility(self,
                          valuation_date,
                          issuer_curve,
                          optionValue):
        """ Calculate the implied CDS option volatility from a price. """
        argtuple = (self, valuation_date, issuer_curve, optionValue)
        sigma = optimize.newton(fvol, x0=0.3, args=argtuple, tol=1e-6,
                                maxiter=50)
        return sigma

##########################################################################
