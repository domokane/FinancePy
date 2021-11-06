###############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
###############################################################################

###############################################################################
# ADD Inflation assumption and resulting yield calculation
# Pricing using inflation curve and discount curve
###############################################################################


from ...utils.date import Date
from ...utils.error import FinError
from ...utils.frequency import annual_frequency, FrequencyTypes
from ...utils.day_count import DayCountTypes
from ...utils.helpers import label_to_string, check_argument_types
from ..bonds.bond import Bond, YTMCalcType

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
                 base_cpi_value: float,
                 num_ex_dividend_days: int = 0):  # Value of CPI index at bond issue date
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
        self._frequency = annual_frequency(freq_type)
        self._face_amount = face_amount  # This is the bond holding size
        self._baseCPIValue = base_cpi_value  # CPI value at issue date of bond
        self._par = 100.0  # This is how price is quoted
        self._redemption = 1.0  # Amount paid at maturity
        self._num_ex_dividend_days = num_ex_dividend_days
        self._inflation_accrued_interest = 0.0

        self._flow_dates = []
        self._flow_amounts = []

        self._settlement_date = Date(1, 1, 1900)
        self._accrued_interest = None
        self._accrued_days = 0.0
        self._alpha = 0.0

        self._calculate_flow_dates()
        self._calculate_flows()

###############################################################################

    def inflation_principal(self,
                            settlement_date: Date,
                            ytm: float,
                            reference_cpi: float,
                            convention: YTMCalcType):
        """ Calculate the principal value of the bond based on the face
        amount and the CPI growth. """

        index_ratio = reference_cpi / self._baseCPIValue
        full_price = self.full_price_from_ytm(settlement_date, ytm, convention)
        principal = full_price * self._face_amount / self._par
        principal = principal - self._accrued_interest
        principal *= index_ratio
        return principal

###############################################################################

    def flat_price_from_yield_to_maturity(self,
                                          settlement_date: Date,
                                          ytm: float,
                                          last_cpn_cpi: float,
                                          convention: YTMCalcType):
        """ Calculate the flat clean price value of the bond based on the clean
        price amount and the CPI growth to the last coupon date. """

        index_ratio = last_cpn_cpi / self._baseCPIValue
        clean_price = self.clean_price_from_ytm(
            settlement_date, ytm, convention)
        flat_price = clean_price * self._face_amount / self._par
        flat_price *= index_ratio
        return flat_price

###############################################################################

    def calc_inflation_accrued_interest(self, settlement_date: Date,
                                        reference_cpi):
        """ Calculate the amount of coupon that has accrued between the
        previous coupon date and the settlement date. This is adjusted by the
        index ratio in line with the CPI growth since the bond base CPI date.
        We assume no ex-dividend period.
        """

        self.calc_accrued_interest(settlement_date)
        index_ratio = reference_cpi/self._baseCPIValue
        self._inflation_accrued_interest = self._accrued_interest * index_ratio
        return self._inflation_accrued_interest

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("ISSUE DATE", self._issue_date)
        s += label_to_string("MATURITY DATE", self._maturity_date)
        s += label_to_string("COUPON", self._coupon)
        s += label_to_string("FREQUENCY", self._freq_type)
        s += label_to_string("ACCRUAL TYPE", self._accrual_type)
        s += label_to_string("FACE AMOUNT", self._face_amount)
        s += label_to_string("BASE CPI VALUE", self._baseCPIValue, "")
        return s

###############################################################################

    def _print(self):
        """ Print a list of the unadjusted coupon payment dates used in
        analytic calculations for the bond. """
        print(self)


###############################################################################
