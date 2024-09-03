##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from math import exp, log, sqrt


from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes, DateGenRuleTypes
from ...utils.day_count import DayCount, DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.global_vars import g_days_in_year
from ...utils.math import ONE_MILLION, INV_ROOT_2_PI, N
from ...utils.error import FinError
from ...products.credit.cds_curve import CDSCurve
from ...products.credit.cds import CDS
from ...utils.helpers import check_argument_types
from ...utils.date import Date
from ...utils.helpers import label_to_string

RPV01_INDEX = 1  # 0 is FULL, 1 is CLEAN

###############################################################################


class CDSIndexOption:
    """Class to manage the pricing and risk management of an option to enter
    into a CDS index. Different pricing algorithms are presented."""

    def __init__(
        self,
        expiry_dt: Date,
        maturity_dt: Date,
        index_cpn: float,
        strike_cpn: float,
        notional: float = ONE_MILLION,
        long_protect: bool = True,
        freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
        dc_type: DayCountTypes = DayCountTypes.ACT_360,
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
    ):
        """Initialisation of the class object. Note that a
        large number of the inputs are set to default values in line with
        the standard contract."""

        check_argument_types(self.__init__, locals())

        if expiry_dt > maturity_dt:
            raise FinError("Expiry date after end date")

        if index_cpn < 0.0:
            raise FinError("Index cpn is negative")

        if strike_cpn < 0.0:
            raise FinError("Index Option strike cpn is negative")

        self.expiry_dt = expiry_dt
        self.maturity_dt = maturity_dt
        self.index_cpn = index_cpn
        self.strike_cpn = strike_cpn
        self.notional = notional
        self.long_protect = long_protect

        self.dc_type = dc_type
        self.dg_type = dg_type
        self.cal_type = cal_type
        self.freq_type = freq_type
        self.bd_type = bd_type

        self.cds_contract = CDS(
            self.expiry_dt,
            self.maturity_dt,
            self.index_cpn,
            1.0,
            self.long_protect,
            self.freq_type,
            self.dc_type,
            self.cal_type,
            self.bd_type,
            self.dg_type,
        )

    ###########################################################################

    def value_adjusted_black(
        self, value_dt, index_curve, index_recovery, libor_curve, sigma
    ):
        """This approach uses two adjustments to black's option pricing
        model to value an option on a CDS index."""

        k = self.strike_cpn
        c = self.index_cpn
        time_to_expiry = (self.expiry_dt - value_dt) / g_days_in_year
        df = libor_curve.df(self.expiry_dt)
        q_expiry_index = index_curve.survival_prob(time_to_expiry)

        cds = CDS(value_dt, self.maturity_dt, k)
        strike_curve = CDSCurve(value_dt, [cds], libor_curve, index_recovery)
        #        qExpiryStrike = strike_curve.surv_prob(time_to_expiry)

        strike_rpv01 = self.cds_contract.risky_pv01(value_dt, strike_curve)[
            "clean_rpv01"
        ]
        index_rpv01 = self.cds_contract.risky_pv01(value_dt, index_curve)[
            "clean_rpv01"
        ]

        s = self.cds_contract.par_spread(value_dt, index_curve)

        fep = df * (1.0 - q_expiry_index) * (1.0 - index_recovery)
        adj_fwd = s + fep / index_rpv01
        adj_strike = c + (k - c) * strike_rpv01 / index_rpv01 / q_expiry_index

        denom = sigma * sqrt(time_to_expiry)
        d1 = log(adj_fwd / adj_strike) + 0.5 * sigma * sigma * time_to_expiry
        d2 = log(adj_fwd / adj_strike) - 0.5 * sigma * sigma * time_to_expiry
        d1 /= denom
        d2 /= denom

        v_pay = (adj_fwd * N(d1) - adj_strike * N(d2)) * index_rpv01
        v_rec = (adj_strike * N(-d2) - adj_fwd * N(-d1)) * index_rpv01

        v_pay *= self.notional
        v_rec *= self.notional

        return (v_pay, v_rec)

    ###########################################################################

    def value_anderson(self, value_dt, issuer_curves, index_recovery, sigma):
        """This function values a CDS index option following approach by
        Anderson (2006). This ensures that a no-arbitrage relationship between
        the constituent CDS contract and the CDS index is enforced. It models
        the forward spread as a log-normally distributed quantity and uses the
        credit triangle to compute the forward RPV01."""

        num_credits = len(issuer_curves)
        time_to_expiry = (self.expiry_dt - value_dt) / g_days_in_year
        #        timeToMaturity = (self.maturity_dt - value_dt) / g_days_in_year
        df_to_expiry = issuer_curves[0].df(time_to_expiry)
        libor_curve = issuer_curves[0].libor_curve

        k = self.strike_cpn
        c = self.index_cpn

        strike_cds = CDS(
            self.expiry_dt, self.maturity_dt, self.strike_cpn, 1.0
        )

        strike_curve = CDSCurve(
            value_dt, [strike_cds], libor_curve, index_recovery
        )
        strike_rpv01s = strike_cds.risky_pv01(value_dt, strike_curve)
        q_to_expiry = strike_curve.survival_prob(time_to_expiry)
        strike_value = (k - c) * strike_rpv01s["clean_rpv01"]
        strike_value /= df_to_expiry * q_to_expiry

        exp_h = 0.0
        h1 = 0.0
        h2 = 0.0

        for i_credit in range(0, num_credits):

            issuer_curve = issuer_curves[i_credit]
            q = issuer_curve.survival_prob(time_to_expiry)
            dh1 = (1.0 - issuer_curve.recovery_rate) * (1.0 - q)

            s = self.cds_contract.par_spread(value_dt, issuer_curve)
            rpv01 = self.cds_contract.risky_pv01(value_dt, issuer_curve)
            dh2 = (s - c) * rpv01["clean_rpv01"] / (df_to_expiry * q_to_expiry)

            h1 = h1 + dh1
            h2 = h2 + dh2

        exp_h = (h1 + h2) / num_credits

        x = self._solve_for_x(
            value_dt, sigma, c, index_recovery, libor_curve, exp_h
        )

        v = self._calc_index_payer_option_price(
            value_dt, x, sigma, c, strike_value, libor_curve, index_recovery
        )

        v = v[1]
        v_pay = v * self.notional
        v_rec = v_pay + (strike_value - exp_h) * df_to_expiry * self.notional
        strike_value *= 10000.0
        x *= 10000.0
        exp_h *= 10000.0
        return v_pay, v_rec, strike_value, x, exp_h

    ###########################################################################

    def _solve_for_x(
        self, value_dt, sigma, index_cpn, index_recovery, libor_curve, exp_h
    ):
        """Function to solve for the arbitrage free"""
        x1 = 0.0
        x2 = 0.9999
        ftol = 1e-8
        j_max = 40
        xacc = 0.000000001
        rtb = 999999

        f = (
            self._calc_obj_func(
                x1, value_dt, sigma, index_cpn, index_recovery, libor_curve
            )
            - exp_h
        )

        fmid = (
            self._calc_obj_func(
                x2, value_dt, sigma, index_cpn, index_recovery, libor_curve
            )
            - exp_h
        )

        if f * fmid >= 0.0:
            raise FinError("Solution not bracketed.")

        if f < 0.0:
            rtb = x1
            dx = x2 - x1
        else:
            rtb = x2
            dx = x1 - x2

        for _ in range(0, j_max):
            dx = dx * 0.5
            xmid = rtb + dx
            fmid = (
                self._calc_obj_func(
                    xmid,
                    value_dt,
                    sigma,
                    index_cpn,
                    index_recovery,
                    libor_curve,
                )
                - exp_h
            )
            if fmid <= 0.0:
                rtb = xmid
            if abs(dx) < xacc or abs(fmid) < ftol:
                return rtb

        return rtb

    ###########################################################################

    def _calc_obj_func(
        self,
        x,
        value_dt,
        sigma,
        index_cpn,  # TODO - do I need this input ?
        index_recovery,
        libor_curve,
    ):
        """An internal function used in the Anderson valuation."""

        # The strike value is not relevant here as we want the zeroth element
        # of the return value
        strike_value = 0.0

        values = self._calc_index_payer_option_price(
            value_dt,
            x,
            sigma,
            self.index_cpn,
            strike_value,
            libor_curve,
            index_recovery,
        )

        return values[0]

    ###########################################################################

    def _calc_index_payer_option_price(
        self,
        value_dt,
        x,
        sigma,
        index_cpn,
        strike_value,
        libor_curve,
        index_recovery,
    ):
        """Calculates the intrinsic value of the index payer swap and the
        value of the index payer option which are both returned in an array.
        """

        z = -6.0
        dz = 0.2
        num_z_steps = int(2.0 * abs(z) / dz)

        flow_dts = self.cds_contract.payment_dts
        num_flows = len(flow_dts)
        t_exp = (self.expiry_dt - value_dt) / g_days_in_year
        df_to_expiry = libor_curve.df(self.expiry_dt)
        lgd = 1.0 - index_recovery

        fwd_dfs = [1.0] * (num_flows)
        expiry_to_flow_times = [1.0] * (num_flows)

        for i_flow in range(0, num_flows):
            expiry_to_flow_times[i_flow] = (
                flow_dts[i_flow] - self.expiry_dt
            ) / g_days_in_year
            fwd_dfs[i_flow] = libor_curve.df(flow_dts[i_flow]) / df_to_expiry

        int_h = 0.0
        int_max_h = 0.0

        day_count = DayCount(self.dc_type)

        #  Previous cpn date is last cpn date before valuation date
        for dt in flow_dts:
            pcd = dt
            if dt > value_dt:
                break

        eff = self.expiry_dt
        accrual_factor_pcd_to_expiry = day_count.year_frac(pcd, eff)[0]

        s0 = exp(-0.5 * sigma * sigma * t_exp)

        for _ in range(0, num_z_steps):
            s = x * s0 * exp(sigma * sqrt(t_exp) * z)
            pdf = exp(-(z**2) / 2.0)
            z = z + dz

            fwd_rpv01 = 0.0
            for i_flow in range(1, num_flows):
                acc_factor = self.cds_contract.accrual_factors[i_flow]
                surv_prob = exp(-s * expiry_to_flow_times[i_flow] / lgd)
                fwd_rpv01 += acc_factor * surv_prob * fwd_dfs[i_flow]

            fwd_rpv01 += -accrual_factor_pcd_to_expiry
            h = (s - index_cpn) * fwd_rpv01
            maxh = max(h - strike_value, 0.0)

            int_h += h * pdf
            int_max_h += maxh * pdf

        int_h *= INV_ROOT_2_PI * dz
        int_max_h *= INV_ROOT_2_PI * dz * df_to_expiry
        return int_h, int_max_h

    ###########################################################################

    def __repr__(self):
        """print out details of the CDS contract and all of the calculated
        cash flows"""
        s = label_to_string("OBJECT TYPE", type(self).__name__)
        s += label_to_string("EXPIRY DATE", self.expiry_dt)
        s += label_to_string("MATURITY DATE", self.maturity_dt)
        s += label_to_string("INDEX cpn", self.index_cpn * 10000, "bp\n")
        s += label_to_string("NOTIONAL", self.notional)
        s += label_to_string("LONG PROTECTION", self.long_protect)
        s += label_to_string("FREQUENCY", self.freq_type)
        s += label_to_string("DAYCOUNT", self.dc_type)
        s += label_to_string("CALENDAR", self.cal_type)
        s += label_to_string("BUSDAYRULE", self.bd_type)
        s += label_to_string("DATEGENRULE", self.dg_type)
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)


###############################################################################
