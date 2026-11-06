##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

from typing import Union

from ....utils import FinError
from ....utils import Date
from ....utils import DayCount, DayCountTypes
from ....utils import FrequencyTypes
from ....utils import CalendarTypes, DateGenRuleTypes
from ....utils import Calendar, BusDayAdjustTypes
from ....utils import Schedule
from ....utils import label_to_string, check_argument_types
from ....utils import ONE_MILLION
from ....utils import SwapTypes
from ....utils.check_values import check_curve_dt

##########################################################################


class IborIborSwap:
    """Class for managing an interest rate basis swap contract. This is a
    contract in which a floating leg with one LIBOR tenor is exchanged for a
    floating leg payment in a different LIBOR tenor. There is no exchange of
    par. The contract is entered into at zero initial cost. The contract lasts
    from a start date to a specified maturity date.

    The value of the contract is the NPV of the two cpn streams. Discounting
    is done on a supplied discount curve which is separate from the discount from
    which the implied index rates are extracted."""

    def __init__(
        self,
        effective_dt: Date,  # Date interest starts to accrue
        term_dt_or_tenor: Union[Date, str],  # Date contract ends
        pay_freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
        pay_dc_type: DayCountTypes = DayCountTypes.THIRTY_E_360,
        rec_freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
        rec_dc_type: DayCountTypes = DayCountTypes.THIRTY_E_360,
        basis_swap_spread: float = 0.0,
        notional: float = ONE_MILLION,
        cal_type: CalendarTypes | list | tuple = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
    ):
        """Create an Ibor basis swap contract giving the contract start
        date, its maturity, frequency and day counts on the two floating
        legs and notional. The floating leg parameters have default
        values that can be overwritten if needed. The start date is contractual
        and is the same as the settlement date for a new swap. It is the date
        on which interest starts to accrue. The end of the contract is the
        termination date. This is not adjusted for business days. The adjusted
        termination date is called the maturity date. This is calculated."""

        check_argument_types(self.__init__, locals())

        if isinstance(term_dt_or_tenor, Date):
            self._termination_dt = term_dt_or_tenor
        else:
            self._termination_dt = effective_dt.add_tenor(term_dt_or_tenor)

        calendar = Calendar(cal_type)
        self.maturity_dt = calendar.adjust(self._termination_dt, bd_type)

        if effective_dt > self.maturity_dt:
            raise FinError("Start date after maturity date")

        self.effective_dt = effective_dt
        self.notional = notional
        self._basis_swap_spread = basis_swap_spread
        self._pay_freq_type = pay_freq_type
        self._rec_freq_type = rec_freq_type
        self._pay_dc_type = pay_dc_type
        self._rec_dc_type = rec_dc_type
        self._cal_type = cal_type
        self._bd_type = bd_type
        self._dg_type = dg_type

        self._pay_float_dts = self._generate_float_leg_dts(self._pay_freq_type)
        self._rec_float_dts = self._generate_float_leg_dts(self._rec_freq_type)

        self._pay_leg = None
        self._rec_leg = None

        self.value_dt = None

    ##########################################################################

    def _generate_float_leg_dts(self, freq_type):
        """Generate the floating leg payment dates all the way back to
        the start date of the swap which may precede the valuation date"""

        float_dts = Schedule(
            self.effective_dt,
            self._termination_dt,
            freq_type,
            self._cal_type,
            self._bd_type,
            self._dg_type,
        ).generate()

        return float_dts

    ##########################################################################

    def value(
        self,
        value_dt,
        discount_curve,
        pay_index_curve,
        rec_index_curve,
        pay_first_fixing_rate=None,
        rec_first_fixing_rate=None,
        principal=0.0,
    ):
        """Value the basis swap on a value date given a discount curve and an
        index curve for each of the two floating legs. The basis swap spread
        is added to the floating rate on the received leg. The value is the
        PV of the received leg less the PV of the paid leg."""

        check_curve_dt(value_dt, discount_curve, pay_index_curve, rec_index_curve)

        self.value_dt = value_dt

        pay_float_leg_value, self._pay_leg = self._float_leg_value(
            value_dt,
            discount_curve,
            pay_index_curve,
            self._pay_float_dts,
            self._pay_dc_type,
            0.0,
            pay_first_fixing_rate,
            principal,
        )

        rec_float_leg_value, self._rec_leg = self._float_leg_value(
            value_dt,
            discount_curve,
            rec_index_curve,
            self._rec_float_dts,
            self._rec_dc_type,
            self._basis_swap_spread,
            rec_first_fixing_rate,
            principal,
        )

        value = rec_float_leg_value - pay_float_leg_value
        return value

    ##########################################################################

    def _float_leg_value(
        self,
        value_dt,  # This should be the settlement date
        discount_curve,
        index_curve,
        float_dts,
        dc_type,
        spread,
        first_fixing_rate,
        principal,
    ):
        """Value one floating leg with payments from an index curve and
        discounting based on a supplied discount curve. The leg is described
        by its payment dates, day count and spread. Returns the leg PV and a
        dictionary of the per-payment details used by the print methods."""

        basis = DayCount(dc_type)

        # The swap may have started in the past but we can only value
        # payments that fall after the valuation date. A payment on the
        # valuation date is treated as already paid, as in SwapFloatLeg.
        start_index = 0
        while start_index < len(float_dts) and float_dts[start_index] <= value_dt:
            start_index += 1

        # If the swap has yet to settle then we do not include the
        # start date of the swap as a cpn payment date.
        if value_dt <= self.effective_dt:
            start_index = 1

        leg = {
            "start_index": start_index,
            "first_fixing_rate": first_fixing_rate,
            "year_fracs": [],
            "flows": [],
            "rates": [],
            "dfs": [],
            "flow_pvs": [],
            "total_pv": [],
        }

        # The swap has matured so there is nothing left to value
        if start_index >= len(float_dts):
            return 0.0, leg

        # All discounting is done forward to the valuation date
        df_value_dt = discount_curve.df(value_dt)

        pv = 0.0
        prev_dt = float_dts[start_index - 1]
        df1_index = None

        for i_flow, next_dt in enumerate(float_dts[start_index:]):

            alpha = basis.year_frac(prev_dt, next_dt)[0]

            # The first floating payment is usually already fixed so is
            # not implied by the index curve. The index curve is only read
            # once a rate has to be implied, so a seasoned swap with a
            # supplied fixing never needs a discount factor from before
            # the valuation date.
            if i_flow == 0 and first_fixing_rate is not None:
                float_rate = first_fixing_rate
            else:
                if df1_index is None:
                    df1_index = index_curve.df(prev_dt)
                df2_index = index_curve.df(next_dt)
                float_rate = (df1_index / df2_index - 1.0) / alpha
                df1_index = df2_index

            flow = (float_rate + spread) * alpha * self.notional
            df_discount = discount_curve.df(next_dt) / df_value_dt
            pv += flow * df_discount

            leg["year_fracs"].append(alpha)
            leg["flows"].append(flow)
            leg["rates"].append(float_rate)
            leg["dfs"].append(df_discount)
            leg["flow_pvs"].append(flow * df_discount)
            leg["total_pv"].append(pv)

            prev_dt = next_dt

        flow = principal * self.notional
        pv = pv + flow * df_discount
        leg["flows"][-1] += flow
        leg["flow_pvs"][-1] += flow * df_discount
        leg["total_pv"][-1] = pv

        return pv, leg

    ##########################################################################

    def print_float_leg_pv(self, leg_type: SwapTypes = SwapTypes.PAY):
        """Prints the dates, accrual factors, discount factors, forward rates,
        implied cash amounts, their present value and their cumulative PV for
        the paid or received floating leg using the last valuation."""

        if leg_type == SwapTypes.PAY:
            leg = self._pay_leg
            float_dts = self._pay_float_dts
            freq_type = self._pay_freq_type
            dc_type = self._pay_dc_type
            spread = 0.0
        else:
            leg = self._rec_leg
            float_dts = self._rec_float_dts
            freq_type = self._rec_freq_type
            dc_type = self._rec_dc_type
            spread = self._basis_swap_spread

        print("START DATE:", self.effective_dt)
        print("MATURITY DATE:", self.maturity_dt)
        print("LEG TYPE:", str(leg_type))
        print("SPREAD cpn (%):", spread * 100)
        print("FLOAT LEG FREQUENCY:", str(freq_type))
        print("FLOAT LEG DAY COUNT:", str(dc_type))
        print("VALUATION DATE", self.value_dt)

        if leg is None or len(leg["flows"]) == 0:
            print("Floating Flows not calculated.")
            return

        if leg["first_fixing_rate"] is None:
            print("         *** FIRST FLOATING RATE PAYMENT IS IMPLIED ***")

        header = "PAYMENT_DT     YEAR_FRAC    RATE(%)       PAYMENT         DF"
        header += "         DF*PAYMENT       CUM_PV"
        print(header)

        # By definition the discount factor is 1.0 on the valuation date
        print("%15s %10s %10s %12s %12.8f %12s %12s" % (self.value_dt, "-", "-", "-", 1.0, "-", "-"))

        start_index = leg["start_index"]

        for i_flow, PAYMENT_DT in enumerate(float_dts[start_index:]):
            print(
                "%15s %10.7f %10.5f %12.2f %12.8f %12.2f %12.2f"
                % (
                    PAYMENT_DT,
                    leg["year_fracs"][i_flow],
                    leg["rates"][i_flow] * 100.0,
                    leg["flows"][i_flow],
                    leg["dfs"][i_flow],
                    leg["flow_pvs"][i_flow],
                    leg["total_pv"][i_flow],
                )
            )

    ##########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("START DATE", self.effective_dt)
        s += label_to_string("TERMINATION DATE", self._termination_dt)
        s += label_to_string("MATURITY DATE", self.maturity_dt)
        s += label_to_string("NOTIONAL", self.notional)
        s += label_to_string("BASIS SWAP SPREAD", self._basis_swap_spread)
        s += label_to_string("PAY FREQUENCY", self._pay_freq_type)
        s += label_to_string("REC FREQUENCY", self._rec_freq_type)
        s += label_to_string("PAY DAY COUNT", self._pay_dc_type)
        s += label_to_string("REC DAY COUNT", self._rec_dc_type)
        s += label_to_string("CALENDAR", self._cal_type)
        s += label_to_string("BUS DAY ADJUST", self._bd_type)
        s += label_to_string("DATE GEN TYPE", self._dg_type)
        return s

    ########################################################################################

    def _print(self):
        """Print a list of the unadjusted cpn payment dates used in
        analytic calculations for the bond."""
        print(self)


########################################################################################
