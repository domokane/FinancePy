##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO: Implied volatility
# TODO: Term structure of volatility
# TODO: Check that curve anchor date is valuation date ?

from typing import Optional

from ...utils.date import Date
from ...utils.calendar import Calendar
from ...utils.calendar import CalendarTypes
from ...utils.calendar import DateGenRuleTypes
from ...utils.calendar import BusDayAdjustTypes
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.global_vars import gDaysInYear
from ...utils.math import ONE_MILLION
from ...utils.error import FinError
from ...utils.schedule import Schedule
from ...utils.helpers import label_to_string, check_argument_types
from ...models.black import Black
from ...models.black_shifted import BlackShifted
from ...models.bachelier import Bachelier
from ...models.sabr import SABR
from ...models.sabr_shifted import SABRShifted
from ...models.hw_tree import HWTree
from ...utils.global_types import FinCapFloorTypes, OptionTypes

##########################################################################

from enum import Enum


class IborCapFloorModelTypes(Enum):
    BLACK = 1
    SHIFTED_BLACK = 2
    SABR = 3

##########################################################################


class IborCapFloor():
    """ Class for Caps and Floors. These are contracts which observe a Ibor
    reset L on a future start date and then make a payoff at the end of the
    Ibor period which is Max[L-K,0] for a cap and Max[K-L,0] for a floor.
    This is then day count adjusted for the Ibor period and then scaled by
    the contract notional to produce a valuation. A number of models can be
    selected from."""

    def __init__(self,
                 start_dt: Date,
                 maturity_dt_or_tenor: (Date, str),
                 option_type: FinCapFloorTypes,
                 strike_rate: float,
                 last_fixing: Optional[float] = None,
                 freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
                 dc_type: DayCountTypes = DayCountTypes.THIRTY_E_360_ISDA,
                 notional: float = ONE_MILLION,
                 cal_type: CalendarTypes = CalendarTypes.WEEKEND,
                 bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
                 dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD):
        """ Initialise IborCapFloor object. """

        check_argument_types(self.__init__, locals())

        self._cal_type = cal_type
        self._bd_type = bd_type

        if type(maturity_dt_or_tenor) == Date:
            maturity_dt = maturity_dt_or_tenor
        else:
            maturity_dt = start_dt.add_tenor(maturity_dt_or_tenor)
            calendar = Calendar(self._cal_type)
            maturity_dt = calendar.adjust(maturity_dt,
                                            self._bd_type)

        if start_dt > maturity_dt:
            raise FinError("Start date must be before maturity date")

        self._start_dt = start_dt
        self._maturity_dt = maturity_dt
        self._option_type = option_type
        self._strike_rate = strike_rate
        self._last_fixing = last_fixing
        self._freq_type = freq_type
        self._dc_type = dc_type
        self._notional = notional
        self._dg_type = dg_type

        self._cap_floor_let_values = []
        self._cap_floor_let_alphas = []
        self._cap_floor_let_fwd_rates = []
        self._cap_floor_let_intrinsic = []
        self._cap_floor_let_dfs = []
        self._cap_floor_pv = []

        self._value_dt = None
        self._day_counter = None

###############################################################################

    def _generate_dts(self):

        schedule = Schedule(self._start_dt,
                            self._maturity_dt,
                            self._freq_type,
                            self._cal_type,
                            self._bd_type,
                            self._dg_type)

        self._capFloorLetDates = schedule._adjusted_dts

##########################################################################

    def value(self, value_dt, libor_curve, model):
        """ Value the cap or floor using the chosen model which specifies
        the volatility of the Ibor rate to the cap start date. """

        self._value_dt = value_dt
        self._generate_dts()

        self._day_counter = DayCount(self._dc_type)
        num_options = len(self._capFloorLetDates)
        strike_rate = self._strike_rate

        if strike_rate < 0.0:
            raise FinError("Strike < 0.0")

        if num_options <= 1:
            raise FinError("Number of options in capfloor equals 1")

        self._cap_floor_let_values = [0]
        self._cap_floor_let_alphas = [0]
        self._cap_floor_let_fwd_rates = [0]
        self._cap_floor_let_intrinsic = [0]
        self._cap_floor_let_dfs = [1.00]
        self._cap_floor_pv = [0.0]

        cap_floor_value = 0.0
        cap_floor_let_value = 0.0
        # Value the first caplet or floorlet with known payoff

        start_dt = self._start_dt
        end_dt = self._capFloorLetDates[1]

        if self._last_fixing is None:
            fwd_rate = libor_curve.fwd_rate(start_dt, end_dt,
                                            self._dc_type)
        else:
            fwd_rate = self._last_fixing

        alpha = self._day_counter.year_frac(start_dt, end_dt)[0]
        df = libor_curve.df(end_dt)

        if self._option_type == FinCapFloorTypes.CAP:
            cap_floor_let_value = df * alpha * max(fwd_rate - strike_rate, 0.0)
        elif self._option_type == FinCapFloorTypes.FLOOR:
            cap_floor_let_value = df * alpha * max(strike_rate - fwd_rate, 0.0)

        cap_floor_let_value *= self._notional
        cap_floor_value += cap_floor_let_value

        self._cap_floor_let_fwd_rates.append(fwd_rate)
        self._cap_floor_let_values.append(cap_floor_let_value)
        self._cap_floor_let_alphas.append(alpha)
        self._cap_floor_let_intrinsic.append(cap_floor_let_value)
        self._cap_floor_let_dfs.append(df)
        self._cap_floor_pv.append(cap_floor_value)

        for i in range(2, num_options):

            start_dt = self._capFloorLetDates[i - 1]
            end_dt = self._capFloorLetDates[i]
            alpha = self._day_counter.year_frac(start_dt, end_dt)[0]

            df = libor_curve.df(end_dt)
            fwd_rate = libor_curve.fwd_rate(start_dt, end_dt,
                                            self._dc_type)

            if self._option_type == FinCapFloorTypes.CAP:
                intrinsic_value = df * alpha * max(fwd_rate - strike_rate, 0.0)
            elif self._option_type == FinCapFloorTypes.FLOOR:
                intrinsic_value = df * alpha * max(strike_rate - fwd_rate, 0.0)

            intrinsic_value *= self._notional

            cap_floor_let_value = self.value_caplet_floor_let(value_dt,
                                                           start_dt,
                                                           end_dt,
                                                           libor_curve,
                                                           model)

            cap_floor_value += cap_floor_let_value

            self._cap_floor_let_fwd_rates.append(fwd_rate)
            self._cap_floor_let_values.append(cap_floor_let_value)
            self._cap_floor_let_alphas.append(alpha)
            self._cap_floor_let_intrinsic.append(intrinsic_value)
            self._cap_floor_let_dfs.append(df)
            self._cap_floor_pv.append(cap_floor_value)

        return cap_floor_value

###############################################################################

    def value_caplet_floor_let(self,
                               value_dt,
                               caplet_start_dt,
                               caplet_end_dt,
                               libor_curve,
                               model):
        """ Value the caplet or floorlet using a specific model. """

        t_exp = (caplet_start_dt - self._start_dt) / gDaysInYear

        alpha = self._day_counter.year_frac(caplet_start_dt,
                                            caplet_end_dt)[0]

        f = libor_curve.fwd_rate(caplet_start_dt, caplet_end_dt,
                                 self._dc_type)

        k = self._strike_rate
        df = libor_curve.df(caplet_end_dt)

        if k == 0.0:
            k = 1e-10

        if isinstance(model, Black):

            if self._option_type == FinCapFloorTypes.CAP:
                cap_floor_let_value = model.value(f, k, t_exp, df,
                                                  OptionTypes.EUROPEAN_CALL)
            elif self._option_type == FinCapFloorTypes.FLOOR:
                cap_floor_let_value = model.value(f, k, t_exp, df,
                                                  OptionTypes.EUROPEAN_PUT)

        elif isinstance(model, BlackShifted):

            if self._option_type == FinCapFloorTypes.CAP:
                cap_floor_let_value = model.value(f, k, t_exp, df,
                                                  OptionTypes.EUROPEAN_CALL)
            elif self._option_type == FinCapFloorTypes.FLOOR:
                cap_floor_let_value = model.value(f, k, t_exp, df,
                                                  OptionTypes.EUROPEAN_PUT)

        elif isinstance(model, Bachelier):

            if self._option_type == FinCapFloorTypes.CAP:
                cap_floor_let_value = model.value(f, k, t_exp, df,
                                                  OptionTypes.EUROPEAN_CALL)
            elif self._option_type == FinCapFloorTypes.FLOOR:
                cap_floor_let_value = model.value(f, k, t_exp, df,
                                                  OptionTypes.EUROPEAN_PUT)

        elif isinstance(model, SABR):

            if self._option_type == FinCapFloorTypes.CAP:
                cap_floor_let_value = model.value(f, k, t_exp, df,
                                                  OptionTypes.EUROPEAN_CALL)
            elif self._option_type == FinCapFloorTypes.FLOOR:
                cap_floor_let_value = model.value(f, k, t_exp, df,
                                                  OptionTypes.EUROPEAN_PUT)

        elif isinstance(model, SABRShifted):

            if self._option_type == FinCapFloorTypes.CAP:
                cap_floor_let_value = model.value(f, k, t_exp, df,
                                                  OptionTypes.EUROPEAN_CALL)
            elif self._option_type == FinCapFloorTypes.FLOOR:
                cap_floor_let_value = model.value(f, k, t_exp, df,
                                                  OptionTypes.EUROPEAN_PUT)

        elif isinstance(model, HWTree):

            t_mat = (caplet_end_dt - value_dt) / gDaysInYear
            alpha = self._day_counter.year_frac(caplet_start_dt,
                                                caplet_end_dt)[0]
            strike_price = 1.0/(1.0 + alpha * self._strike_rate)
            notional_adj = (1.0 + self._strike_rate * alpha)
            face_amount = 1.0
            df_times = libor_curve._times
            df_values = libor_curve._dfs

            v = model.option_on_zcb(t_exp, t_mat, strike_price, face_amount,
                                    df_times, df_values)

            # we divide by alpha to offset the multiplication above
            if self._option_type == FinCapFloorTypes.CAP:
                cap_floor_let_value = v['put'] * notional_adj / alpha
            elif self._option_type == FinCapFloorTypes.FLOOR:
                cap_floor_let_value = v['call'] * notional_adj / alpha

        else:
            raise FinError("Unknown model type " + str(model))

        cap_floor_let_value *= (self._notional * alpha)

        return cap_floor_let_value

###############################################################################

    def print_leg(self):
        """ Prints the cap floor payment amounts. """

        print("START DATE:", self._start_dt)
        print("MATURITY DATE:", self._maturity_dt)
        print("OPTION TYPE", str(self._option_type))
        print("STRIKE (%):", self._strike_rate * 100)
        print("FREQUENCY:", str(self._freq_type))
        print("DAY COUNT:", str(self._dc_type))
        print("VALUATION DATE", self._value_dt)

        if len(self._cap_floor_let_values) == 0:
            print("Caplets not calculated.")
            return

        if self._option_type == FinCapFloorTypes.CAP:
            header = "PAYMENT_dt     YEAR_FRAC   FWD_RATE    INTRINSIC      "
            header += "     DF    CAPLET_PV       CUM_PV"
        elif self._option_type == FinCapFloorTypes.FLOOR:
            header = "PAYMENT_dt     YEAR_FRAC   FWD_RATE    INTRINSIC      "
            header += "     DF    FLRLET_PV       CUM_PV"

        print(header)

        iFlow = 0

        for payment_dt in self._capFloorLetDates[iFlow:]:
            if iFlow == 0:
                print("%15s %10s %9s %12s %12.6f %12s %12s" %
                      (payment_dt, "-", "-", "-",
                       self._cap_floor_let_dfs[iFlow], "-", "-"))
            else:
                print("%15s %10.7f %9.5f %12.2f %12.6f %12.2f %12.2f" %
                      (payment_dt,
                       self._cap_floor_let_alphas[iFlow],
                       self._cap_floor_let_fwd_rates[iFlow]*100,
                       self._cap_floor_let_intrinsic[iFlow],
                       self._cap_floor_let_dfs[iFlow],
                       self._cap_floor_let_values[iFlow],
                       self._cap_floor_pv[iFlow]))

            iFlow += 1

###############################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("START DATE", self._start_dt)
        s += label_to_string("MATURITY DATE", self._maturity_dt)
        s += label_to_string("STRIKE COUPON", self._strike_rate * 100)
        s += label_to_string("OPTION TYPE", str(self._option_type))
        s += label_to_string("FREQUENCY", str(self._freq_type))
        s += label_to_string("DAY COUNT", str(self._dc_type), "")
        return s

###############################################################################

    def _print(self):
        print(self)

###############################################################################
