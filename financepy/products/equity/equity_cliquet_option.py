##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from ...utils.frequency import FrequencyTypes
from ...utils.error import FinError
from ...utils.global_types import OptionTypes

from ...utils.helpers import label_to_string, check_argument_types
from ...utils.helpers import option_years
from ...utils.date import Date
from ...utils.calendar import BusDayAdjustTypes
from ...utils.calendar import CalendarTypes, DateGenRuleTypes
from ...utils.schedule import Schedule
from ...utils.check_values import check_curve_dt

from ...products.equity.equity_option import EquityOption
from ...market.curves.flat_discount_curve import DiscountCurve

from ...models.black_scholes_analytic import european_value
from ...models.black_scholes import BlackScholes
from ...models.model import Model

########################################################################################
# TODO: Do we need to day count adjust option payoffs ?
# TODO: Monte Carlo pricer
########################################################################################


class EquityCliquetOption(EquityOption):
    """A EquityCliquetOption is a series of options which start and stop at
    successive times with each subsequent option resetting its strike to be ATM
    at the start of its life. This is also known as a reset option."""

    def __init__(
        self,
        start_dt: Date,
        final_expiry_dt: Date,
        opt_type: OptionTypes,
        freq_type: FrequencyTypes,
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
    ):
        """Create the EquityCliquetOption by passing in the start date
        and the end date and whether it is a call or a put. Some additional
        data is needed in order to calculate the individual payments."""

        check_argument_types(self.__init__, locals())

        if opt_type != OptionTypes.EUROPEAN_CALL and opt_type != OptionTypes.EUROPEAN_PUT:
            raise FinError("Unknown OPTION_TYPE" + str(opt_type))

        if final_expiry_dt < start_dt:
            raise FinError("Expiry date precedes start date")

        self.start_dt = start_dt
        self.final_expiry_dt = final_expiry_dt
        self.opt_type = opt_type
        self.freq_type = freq_type
        self.cal_type = cal_type
        self.bd_type = bd_type
        self.dg_type = dg_type

        self.v_options = None
        self._dfs = None
        self.actual_dts = None

        self.expiry_dts = Schedule(
            self.start_dt,
            self.final_expiry_dt,
            self.freq_type,
            self.cal_type,
            self.bd_type,
            self.dg_type,
        ).generate()

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
    ):
        """Value the cliquet option as a sequence of options using the Black-
        Scholes model."""

        if isinstance(value_dt, Date) is False:
            raise FinError("Valuation date is not a Date")

        if value_dt > self.final_expiry_dt:
            raise FinError("Value date after final expiry date.")

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        s0 = stock_price
        v_cliquet = 0.0

        self.v_options = []
        self.dfs = []
        self.actual_dts = []

        if isinstance(model, BlackScholes):

            fwd_vol = model.volatility
            fwd_vol = max(fwd_vol, 1e-6)

            dt_prev = value_dt

            for dt in self.expiry_dts:

                if dt > value_dt:

                    t_vol = option_years(dt_prev, dt)

                    df_end = discount_curve.df(dt)
                    dq_start = dividend_curve.df(dt_prev)

                    # The deflator is out to the option reset time
                    fwd_r = discount_curve.fwd_zero_rate_cc(dt_prev, dt)
                    fwd_q = dividend_curve.fwd_zero_rate_cc(dt_prev, dt)

                    v = european_value(1.0, t_vol, 1.0, fwd_r, fwd_q, fwd_vol, self.opt_type.value)
                    v_fwd_opt = s0 * dq_start * v
                    v_cliquet += v_fwd_opt

                    self.dfs.append(df_end)
                    self.v_options.append(v_fwd_opt)
                    self.actual_dts.append(dt)

                    dt_prev = dt
        else:
            raise FinError("Unknown Model Type")

        return v_cliquet

    ###########################################################################

    def print_payments(self):

        if self.v_options is None:
            raise FinError("Options not created yet.")

        num_options = len(self.v_options)
        for i in range(0, num_options):
            print(self.actual_dts[i], self.dfs[i], self.v_options[i])

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("START_DATE", self.start_dt)
        s += label_to_string("FINAL EXPIRY DATE", self.final_expiry_dt)
        s += label_to_string("OPTION_TYPE", self.opt_type)
        s += label_to_string("FREQUENCY TYPE", self.freq_type)
        s += label_to_string("CALENDAR TYPE", self.cal_type)
        s += label_to_string("BUS_DAY_ADJUST", self.bd_type)
        s += label_to_string("DATE GEN RULE TYPE", self.dg_type, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
