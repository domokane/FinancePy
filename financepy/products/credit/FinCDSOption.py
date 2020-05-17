# -*- coding: utf-8 -*-
"""
Created on Sat Aug  3 14:15:33 2019

@author: Dominic
"""
from math import sqrt, log
from scipy import optimize

from ...finutils.FinCalendar import FinCalendarTypes
from ...finutils.FinCalendar import FinDayAdjustTypes, FinDateGenRuleTypes
from ...finutils.FinDayCount import FinDayCountTypes
from ...finutils.FinFrequency import FinFrequencyTypes
from ...finutils.FinGlobalVariables import gDaysInYear
from ...finutils.FinMath import ONE_MILLION, N
from ...products.credit.FinCDS import FinCDS

##########################################################################


def fvol(volatility, *args):
    ''' Root searching function in the calculation of the CDS implied volatility. '''

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
    ''' Class to manage the pricing and risk-management of options on a CDS. '''

    def __init__(self,
                 expiryDate,
                 maturityDate,
                 strikeCoupon,
                 notional=ONE_MILLION,
                 longProtection=True,
                 knockoutFlag=True,
                 frequencyType=FinFrequencyTypes.QUARTERLY,
                 dayCountType=FinDayCountTypes.ACT_360,
                 calendarType=FinCalendarTypes.WEEKEND,
                 busDayAdjustType=FinDayAdjustTypes.FOLLOWING,
                 dateGenRuleType=FinDateGenRuleTypes.BACKWARD):

        if maturityDate < expiryDate:
            raise ValueError("Maturity date must be after option expiry date")

        if strikeCoupon < 0.0:
            raise ValueError("Strike must be greater than zero")

        if frequencyType not in FinFrequencyTypes:
            raise ValueError(
                "Unknown Fixed Frequency type " +
                str(frequencyType))

        if calendarType not in FinCalendarTypes:
            raise ValueError("Unknown Calendar type " + str(calendarType))

        if busDayAdjustType not in FinDayAdjustTypes:
            raise ValueError(
                "Unknown Business Day Adjust type " +
                str(busDayAdjustType))

        if dateGenRuleType not in FinDateGenRuleTypes:
            raise ValueError(
                "Unknown Date Gen Rule type " +
                str(dateGenRuleType))

        self._expiryDate = expiryDate
        self._maturityDate = maturityDate
        self._strikeCoupon = strikeCoupon
        self._longProtection = longProtection
        self._knockoutFlag = knockoutFlag
        self._notional = notional

        self._frequencyType = frequencyType
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
            raise ValueError("Expiry date is now or in the past")

        if volatility < 0.0:
            raise ValueError("Volatility must be greater than zero")

        # The underlying is a forward starting option that steps in on
        # the expiry date and matures on the expiry date with a coupon
        # set equal to the option spread strike
        cds = FinCDS(self._expiryDate,
                     self._maturityDate,
                     self._strikeCoupon,
                     self._notional,
                     self._longProtection,
                     self._frequencyType,
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
        if self._knockoutFlag == False and self._longProtection == True:

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
