###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

###############################################################################
# ADD Inflation assumption and resulting yield calculation
# Pricing using inflation curve and discount curve
###############################################################################


from ...utils.date import Date
from ...utils.FinError import FinError
from ...utils.frequency import Frequency, FrequencyTypes
from ...utils.day_count import DayCountTypes
from ...utils.helper_functions import labelToString, check_argument_types
from ..bonds.bond import Bond, FinYTMCalcType

###############################################################################


class FinInflationBond(Bond):
    """ Class for inflation-linked bonds like TIPS and related analytics. These
    are bonds with coupon and principal adjusted by an index such as the CPI.
    We inherit from the Bond class. """

    def __init__(self,
                 issue_date: Date,
                 maturity_date: Date,
                 coupon: float,  # Annualised bond coupon before inflation
                 freq_type: FrequencyTypes,
                 accrual_type: DayCountTypes,
                 face_amount: float,
                 baseCPIValue: float,
                 numExDividendDays: int = 0): # Value of CPI index at bond issue date
        """ Create FinInflationBond object by providing Maturity, Frequency,
        coupon, frequency and the accrual convention type. You must also supply
        the base CPI used for all coupon and principal related calculations. 
        The class inherits from Bond so has many similar functions. The YTM"""
        
        check_argument_types(self.__init__, locals())

        if issue_date >= maturity_date:
            raise FinError("Issue Date must preceded maturity date.")

        self._issue_date = issue_date
        self._maturity_date = maturity_date
        self._coupon = coupon
        self._freq_type = freq_type
        self._accrual_type = accrual_type
        self._frequency = Frequency(freq_type)
        self._face_amount = face_amount  # This is the bond holding size
        self._baseCPIValue = baseCPIValue # CPI value at issue date of bond
        self._par = 100.0  # This is how price is quoted
        self._redemption = 1.0 # Amount paid at maturity
        self._numExDividendDays = numExDividendDays

        self._flow_dates = []
        self._flow_amounts = []

        self._settlement_date = Date(1, 1, 1900)
        self._accruedInterest = None
        self._accrued_days = 0.0
        self._alpha = 0.0
                   
        self._calculate_flow_dates()
        self._calculateFlows()

###############################################################################

    def inflationPrincipal(self,
                           settlement_date: Date,
                           ytm: float,
                           referenceCPI: float,
                           convention: FinYTMCalcType):
        """ Calculate the principal value of the bond based on the face
        amount and the CPI growth. """

        indexRatio = referenceCPI / self._baseCPIValue
        full_price = self.full_price_from_ytm(settlement_date, ytm, convention)
        principal = full_price * self._face_amount / self._par
        principal = principal - self._accruedInterest
        principal *= indexRatio
        return principal

###############################################################################

    def flatPriceFromYieldToMaturity(self,
                                     settlement_date: Date,
                                     ytm: float,
                                     lastCpnCPI: float,
                                     convention: FinYTMCalcType):
        """ Calculate the flat clean price value of the bond based on the clean
        price amount and the CPI growth to the last coupon date. """

        indexRatio = lastCpnCPI / self._baseCPIValue
        clean_price = self.clean_price_from_ytm(settlement_date, ytm, convention)
        flatPrice = clean_price * self._face_amount / self._par
        flatPrice *= indexRatio
        return flatPrice

###############################################################################

    def calcInflationAccruedInterest(self, settlement_date: Date,
                                     referenceCPI):
        """ Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date. This is adjusted by the
        index ratio in line with the CPI growth since the bond base CPI date.
        We assume no ex-dividend period.
        """

        self.calc_accrued_interest(settlement_date)
        indexRatio = referenceCPI/self._baseCPIValue
        self._inflationAccruedInterest = self._accruedInterest * indexRatio        
        return self._inflationAccruedInterest

###############################################################################

    def __repr__(self):
        s = labelToString("OBJECT TYPE", type(self).__name__)
        s += labelToString("ISSUE DATE", self._issue_date)
        s += labelToString("MATURITY DATE", self._maturity_date)
        s += labelToString("COUPON", self._coupon)
        s += labelToString("FREQUENCY", self._freq_type)
        s += labelToString("ACCRUAL TYPE", self._accrual_type)
        s += labelToString("FACE AMOUNT", self._face_amount)
        s += labelToString("BASE CPI VALUE", self._baseCPIValue, "")
        return s

###############################################################################

    def _print(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        print(self)


###############################################################################
