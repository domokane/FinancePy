##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from math import sqrt, log
from scipy import optimize


from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinBusDayAdjustTypes, FinDateGenRuleTypes
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinMath import ONE_MILLION, N
from ...products.credit.FinCDS import FinCDS
from ...finutils.FinHelperFunctions import checkArgumentTypes 
from ...finutils.FinDate import FinDate
from ...finutils.FinError import FinError

##########################################################################


def fvol(volatility, *args):
    ''' Root searching function in the calculation of the CDS implied
    volatility. '''

    self = args[0]
    valuationDate = args[1]
    issuerCurve = args[2]
    optionPrice = args[3]
    value = self.value(valuationDate, issuerCurve, volatility)
    objFn = value - optionPrice

    return objFn

##########################################################################
##########################################################################


class FinCDSOption():
    ''' Class to manage the pricing and risk-management of an option on a
    single-name CDS. This is a contract in which the option buyer pays for an
    option to either buy or sell protection on the underlying CDS at a fixed
    spread agreed today and to be exercised in the future on a specified expiry
    date. The option may or may not cancel if there is a credit event before
    option expiry. This needs to be specified. '''

    def __init__(self,
                 expiryDate: FinDate,
                 maturityDate: FinDate,
                 strikeCoupon: float,
                 notional: float = ONE_MILLION,
                 longProtection: bool = True,
                 knockoutFlag: bool = True,
                 freqType: FinFrequencyTypes = FinFrequencyTypes.QUARTERLY,
                 dayCountType: FinDayCountTypes = FinDayCountTypes.ACT_360,
                 calendarType: FinCalendarTypes = FinCalendarTypes.WEEKEND,
                 busDayAdjustType: FinBusDayAdjustTypes = FinBusDayAdjustTypes.FOLLOWING,
                 dateGenRuleType: FinDateGenRuleTypes = FinDateGenRuleTypes.BACKWARD):
        ''' Create a FinCDSOption object with the option expiry date, the
        maturity date of the underlying CDS, the option strike coupon,
        notional, whether the option knocks out or not in the event of a credit
        event before expiry and the payment details of the underlying CDS. '''

        checkArgumentTypes(self.__init__, locals())

        if maturityDate < expiryDate:
            raise FinError("Maturity date must be after option expiry date")

        if strikeCoupon < 0.0:
            raise FinError("Strike must be greater than zero")

        self._expiryDate = expiryDate
        self._maturityDate = maturityDate
        self._strikeCoupon = strikeCoupon
        self._longProtection = longProtection
        self._knockoutFlag = knockoutFlag
        self._notional = notional

        self._freqType = freqType
        self._dayCountType = dayCountType
        self._calendarType = calendarType
        self._businessDateAdjustType = busDayAdjustType
        self._dateGenRuleType = dateGenRuleType

##########################################################################

    def value(self,
              valuationDate,
              issuerCurve,
              volatility):
        ''' Value the CDS option using Black's model with an adjustment for any
        Front End Protection.
        TODO - Should the CDS be created in the init method ? '''

        if valuationDate > self._expiryDate:
            raise FinError("Expiry date is now or in the past")

        if volatility < 0.0:
            raise FinError("Volatility must be greater than zero")

        # The underlying is a forward starting option that steps in on
        # the expiry date and matures on the expiry date with a coupon
        # set equal to the option spread strike
        cds = FinCDS(self._expiryDate,
                     self._maturityDate,
                     self._strikeCoupon,
                     self._notional,
                     self._longProtection,
                     self._freqType,
                     self._dayCountType,
                     self._calendarType,
                     self._businessDateAdjustType,
                     self._dateGenRuleType)

        strike = self._strikeCoupon
        forwardSpread = cds.parSpread(valuationDate, issuerCurve)
        forwardRPV01 = cds.riskyPV01(valuationDate, issuerCurve)['full_rpv01']

        timeToExpiry = (self._expiryDate - valuationDate) / gDaysInYear
        logMoneyness = log(forwardSpread / strike)

        halfVolSquaredT = 0.5 * volatility * volatility * timeToExpiry
        volSqrtT = volatility * sqrt(timeToExpiry)

        d1 = (logMoneyness + halfVolSquaredT) / volSqrtT
        d2 = (logMoneyness - halfVolSquaredT) / volSqrtT

        if self._longProtection:
            optionValue = forwardSpread * N(d1) - strike * N(d2)
        else:
            optionValue = strike * N(-d2) - forwardSpread * N(-d1)

        optionValue = optionValue * forwardRPV01

        # If the option does not knockout on a default before expiry then we
        # need to include the cost of protection which is provided between
        # the value date and the expiry date
        if self._knockoutFlag is False and self._longProtection is True:

            df = issuerCurve.getDF(timeToExpiry)
            q = issuerCurve.getSurvProb(timeToExpiry)
            recovery = issuerCurve._recoveryRate
            frontEndProtection = df * (1.0 - q) * (1.0 - recovery)
            optionValue += frontEndProtection

        # we return the option price in dollars
        return optionValue * self._notional

##########################################################################

    def impliedVolatility(self,
                          valuationDate,
                          issuerCurve,
                          optionValue):
        ''' Calculate the implied CDS option volatility from a price. '''
        argtuple = (self, valuationDate, issuerCurve, optionValue)
        sigma = optimize.newton(fvol, x0=0.3, args=argtuple, tol=1e-6,
                                maxiter=50)
        return sigma

##########################################################################
