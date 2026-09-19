##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


import numpy as np

from ...utils.frequency import FrequencyTypes
from ...utils.error import FinError
from ...utils.global_types import OptionTypes

from ...utils.helpers import label_to_string, check_argument_types
from ...utils.date import Date
from ...utils.schedule import Schedule
from ...products.equity.equity_option import EquityOption
from ...market.curves.flat_discount_curve import DiscountCurve
from ...utils.check_values import check_curve_dt
from ...utils.helpers import option_years

from ...models.black_scholes_analytic import european_value
from ...models.black_scholes import BlackScholes
from ...models.model import Model

########################################################################################
# TODO: Do we need to day count adjust option payoffs ?
# TODO: Monte Carlo pricer
########################################################################################


class EquityForwardStartOption(EquityOption):
    """A EquityCliquetOption is an option which starts on t_F and sets the strike to
    S(T_F). It expires at time T. This is also known as a reset option."""

    def __init__(self, start_dt: Date, final_expiry_dt: Date, opt_type: OptionTypes, freq_type: FrequencyTypes):
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

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        if value_dt > self.final_expiry_dt:
            raise FinError("Value date after final expiry date.")

        s = stock_price
        v_cliquet = 0.0

        self.v_options = []
        self._dfs = []
        self.actual_dts = []

        call_type = OptionTypes.EUROPEAN_CALL
        put_type = OptionTypes.EUROPEAN_PUT

        if isinstance(model, BlackScholes):

            v = model.volatility
            v = max(v, 1e-6)
            t_prev = 0.0
            dt_prev = None

            for dt in self.expiry_dts:

                if dt > value_dt:

                    df = discount_curve.df(dt)
                    r = discount_curve.zero_rate_cc(dt)
                    t_exp = option_years(value_dt, dt)

                    # option life
                    tau = t_exp - t_prev

                    # The deflator is out to the option reset time
                    dq = dividend_curve.df(dt_prev)

                    # The option dividend is over the option life
                    dq_mat = dividend_curve.df(dt)

                    old_q = -np.log(dq_mat / dq) / tau
                    q = dividend_curve.fwd_zero_rate_cc(dt_prev, dt)
                    print("q old", old_q, q)

                    if self.opt_type == call_type:
                        v_call = european_value(1.0, tau, 1.0, r, q, v, call_type.value)
                        v_fwd_opt = s * dq * v_call
                        v_cliquet += v_fwd_opt
                    elif self.opt_type == put_type:
                        v_put = european_value(1.0, tau, 1.0, r, q, v, put_type.value)
                        v_fwd_opt = s * dq * v_put
                        v_cliquet += v_fwd_opt
                    else:
                        raise FinError("Unknown OPTION_TYPE")

                    #  print(dt, r, df, q, v_fwd_opt, v_cliquet)

                    self._dfs.append(df)
                    self.v_options.append(v)
                    self.actual_dts.append(dt)
                    t_prev = t_exp
                    dt_prev = dt
        else:
            raise FinError("Unknown Model Type")

        return v_cliquet

    ###########################################################################

    def print_payments(self):
        num_options = len(self.v_options)
        for i in range(0, num_options):
            print(self.actual_dts[i], self._dfs[i], self.v_options[i])

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("START_DATE", self.start_dt)
        s += label_to_string("FINAL EXPIRY DATE", self.final_expiry_dt)
        s += label_to_string("OPTION_TYPE", self.opt_type)
        s += label_to_string("FREQUENCY TYPE", self.freq_type)
        s += label_to_string("DC_TYPE", self.accrual_dc_type)
        s += label_to_string("CALENDAR TYPE", self.cal_type)
        s += label_to_string("BUS_DAY_ADJUST", self.bd_type)
        s += label_to_string("DATE GEN RULE TYPE", self.dg_type, "")
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


########################################################################################
