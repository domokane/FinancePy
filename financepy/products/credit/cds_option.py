##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from math import sqrt, log
from scipy import optimize

from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes, DateGenRuleTypes
from ...utils.day_count import DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.global_vars import ONE_MILLION, TO_BPS
from ...utils.math import normcdf
from ...products.credit.cds import CDS
from ...utils.helpers import check_argument_types
from ...utils.date import Date
from ...utils.error import FinError
from ...utils.check_values import check_curve_dt
from ...market.curves.discount_curve import DiscountCurve
from ...market.curves.cds_curve import CDSCurve
from ...utils.helpers import option_years
from ...utils.valuation_log import log_value, log_table, log_row, log_section, log_bool
from ...utils.valuation_log import log_amount, log_float, log_percentage, log_date, log_bps

from ...utils.helpers import label_to_string, table_to_string

########################################################################################


def _fvol(volatility: float, *args):
    """Root searching function in the calculation of the CDS implied volatility."""

    self = args[0]
    value_dt = args[1]
    issuer_curve = args[2]
    option_value = args[3]
    value = self.value(value_dt, issuer_curve, volatility)
    obj_fn = value - option_value
    return obj_fn


########################################################################################


class CDSOption:
    """Class to manage the pricing and risk-management of an option on a
    single-name CDS. This is a contract in which the option buyer pays for an
    option to either buy or sell protection on the underlying CDS at a fixed
    spread agreed today and to be exercised in the future on a specified expiry
    date. The option may or may not knockout if there is a credit event before
    option expiry. This needs to be specified. Default is to knockout on default. """

    def __init__(
        self,
        expiry_dt: Date,
        maturity_dt: Date,
        strike_cpn: float,
        notional: float = ONE_MILLION,
        is_payer: bool = True,
        knockout_flag: bool = True,
        freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
        accrual_dc_type: DayCountTypes = DayCountTypes.ACT_360,
        cal_type: CalendarTypes | list | tuple = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
    ) -> None:
        """Create a CDSOption object with the option expiry date, the
        maturity date of the underlying CDS, the option strike coupon,
        notional, whether the option knocks out or not in the event of a credit
        event before expiry and the payment details of the underlying CDS."""

        check_argument_types(self.__init__, locals())

        if maturity_dt < expiry_dt:
            raise FinError("Maturity date must be after option expiry date")

        if strike_cpn < 0.0:
            raise FinError("Strike must be greater than zero")

        self.expiry_dt = expiry_dt
        self.maturity_dt = maturity_dt
        self.strike_cpn = strike_cpn
        self.is_payer = is_payer
        self.knockout_flag = knockout_flag
        self.notional = notional

        self.freq_type = freq_type
        self.accrual_dc_type = accrual_dc_type
        self.cal_type = cal_type
        self.bd_type = bd_type
        self.dg_type = dg_type

        # The underlying is a forward starting option that steps in on
        # the expiry date and matures on the expiry date with a coupon
        # set equal to the option spread strike
        self.cds = CDS(
            self.expiry_dt,
            self.maturity_dt,
            self.strike_cpn,
            self.notional,
            self.is_payer,  # Same as long protection
            self.freq_type,
            self.accrual_dc_type,
            self.cal_type,
            self.bd_type,
            self.dg_type,
        )

    ####################################################################################

    def value(self,
              value_dt: Date,
              issuer_curve: CDSCurve,
              volatility: float):
        """Value the CDS option using Black's model with an adjustment
        for any Front End Protection."""

        if value_dt > self.expiry_dt:
            raise FinError("Expiry date is now or in the past")

        if volatility < 0.0:
            raise FinError("Volatility must be greater than zero")

        strike = self.strike_cpn

        log_date("Value Date", value_dt, level=2)
        log_date("Expiry Date", self.expiry_dt, level=2)
        log_date("Maturity Date", self.maturity_dt, level=2)
        log_bps("Strike", strike, level=2)
        log_bool("Payer", self.is_payer, level=2)
        log_bool("Knockout", self.knockout_flag, level=2)
        log_amount("Notional", self.notional, level=1)

        check_curve_dt(value_dt, issuer_curve)

        log_section("Valuation Steps", level=2)

        log_percentage("Volatility", volatility, level=2)

        fwd_spread, fwd_rpv01, _ = self.cds.par_spread(value_dt,
                                                       issuer_curve,
                                                       all_values=True)

        t_exp = option_years(value_dt, self.expiry_dt)

        log_bps("Forward Spread", fwd_spread, level=2)
        log_float("Forward RPV01", fwd_rpv01, level=2)
        log_float("Time to Expiry", t_exp, level=2)

        half_vol_squared_t = 0.5 * volatility * volatility * t_exp
        vol_sqrt_t = volatility * sqrt(t_exp)

        # Black formula's zero-variance limit
        if vol_sqrt_t == 0.0:
            if self.long_protection:
                option_value = max(fwd_spread - strike, 0.0)
            else:
                option_value = max(strike - fwd_spread, 0.0)

        else:
            log_moneyness = log(fwd_spread / strike)
            d1 = (log_moneyness + half_vol_squared_t) / vol_sqrt_t
            d2 = (log_moneyness - half_vol_squared_t) / vol_sqrt_t

            log_float("Log Moneyness", log_moneyness)
            log_float("d1", d1)
            log_float("d2", d2)

            if self.is_payer:
                option_value = (
                    fwd_spread * normcdf(d1)
                    - strike * normcdf(d2)
                )
            else:
                option_value = (
                    strike * normcdf(-d2)
                    - fwd_spread * normcdf(-d1)
                )

        option_value = option_value * fwd_rpv01

        log_float("Value Per Unit Notional", option_value, level=1)

        # Front-end protection
        front_end_protection = 0.0

        if self.knockout_flag is False and self.is_payer is True:

            df = issuer_curve.libor_curve.df_t(t_exp)
            q = issuer_curve.survival_prob(t_exp)
            recovery = issuer_curve.recovery_rate

            front_end_protection = df * (1.0 - q) * (1.0 - recovery)

            log_float("FEP Discount Factor", df, level=2)
            log_float("FEP Survival Probability", q, level=2)
            log_percentage("Recovery Rate", recovery, level=2)
            log_float("FEP Per Unit Notional", front_end_protection)

            option_value += front_end_protection

        option_value = option_value * self.notional

        log_amount("Final Option Value", option_value, level=1)

        option_value_in_cts = option_value / self.notional * 10000

        log_amount("Final Option Value (cts)", option_value_in_cts, level=1)

        return option_value

    ####################################################################################

    def implied_volatility(self, value_dt: Date, issuer_curve: DiscountCurve, option_value):
        """Calculate the implied CDS option volatility from a price."""
        arg_tuple = (self, value_dt, issuer_curve, option_value)
        sigma = optimize.newton(_fvol, x0=0.3, args=arg_tuple, tol=1e-6, maxiter=50)
        return sigma

    ###########################################################################

    def __repr__(self):
        """print out details of the CDS contract and all the calculated
        cash flows"""
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EFFECTIVE DATE", self.step_in_dt)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("TERMINATION DATE", self.maturity_dt)
        s += label_to_string("STRIKE SPREAD", self.strike_cpn * 10000, "bp\n")
        s += label_to_string("KNOCKOUT FLAG", self.knockout_flag)
        s += label_to_string("PAYER", self.is_payer)
        s += label_to_string("NOTIONAL", self.notional)
        s += label_to_string("DAY COUNT", self.accrual_dc_type)
        s += label_to_string("FREQUENCY", self.freq_type)
        s += label_to_string("CALENDAR", self.cal_type)
        s += label_to_string("BUS DAY RULE", self.bd_type)
        s += label_to_string("DATE GEN RULE", self.dg_type)

        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################

########################################################################################
