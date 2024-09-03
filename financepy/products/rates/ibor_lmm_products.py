##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################

# TODO: Extend to allow term structure of volatility
# TODO: Extend to allow two fixed legs in underlying swap
# TODO: Cash settled swaptions

""" This module implements the LMM in the spot measure. It combines both model
and product specific code - I am not sure if it is better to separate these. At
the moment this seems to work ok.

THIS IS STILL IN PROTOPTYPE MODE. DO NOT USE. """

import numpy as np

from ...utils.calendar import CalendarTypes
from ...utils.calendar import BusDayAdjustTypes
from ...utils.calendar import DateGenRuleTypes
from ...utils.day_count import DayCountTypes
from ...utils.frequency import FrequencyTypes
from ...utils.day_count import DayCount
from ...utils.schedule import Schedule
from ...utils.error import FinError
from ...utils.helpers import check_argument_types
from ...utils.date import Date

from ...models.lmm_mc import lmm_simulate_fwds_1f
from ...models.lmm_mc import lmm_simulate_fwds_mf
from ...models.lmm_mc import lmm_simulate_fwds_nf
from ...models.lmm_mc import ModelLMMModelTypes
from ...models.lmm_mc import lmm_cap_flr_pricer

from ...utils.global_vars import g_days_in_year
from ...utils.math import ONE_MILLION

from ...utils.global_types import SwapTypes
from ...utils.global_types import FinCapFloorTypes

from financepy.market.volatility.ibor_cap_vol_curve import IborCapVolCurve

###############################################################################


class IborLMMProducts:
    """This is the class for pricing Ibor products using the LMM."""

    def __init__(
        self,
        settle_dt: Date,
        maturity_dt: Date,
        float_freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
        float_dc_type: DayCountTypes = DayCountTypes.THIRTY_E_360,
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
    ):
        """Create a European-style swaption by defining the exercise date of
        the swaption, and all of the details of the underlying interest rate
        swap including the fixed cpn and the details of the fixed and the
        floating leg payment schedules."""

        check_argument_types(self.__init__, locals())

        if settle_dt > maturity_dt:
            raise FinError("Settlement date must be before maturity date")

        """ Set up the grid for the Ibor rates that are to be simulated. These
        must be consistent with the floating rate leg of the product that is to
        be priced. """

        self.start_dt = settle_dt
        self.grid_dts = Schedule(
            settle_dt, maturity_dt, float_freq_type, cal_type, bd_type, dg_type
        ).generate()

        self.accrual_factors = []
        self.float_dc_type = float_dc_type

        basis = DayCount(self.float_dc_type)
        prev_dt = self.grid_dts[0]

        self.grid_times = [0.0]

        for next_dt in self.grid_dts[1:]:
            tau = basis.year_frac(prev_dt, next_dt)[0]
            t = (next_dt - self.grid_dts[0]) / g_days_in_year
            self.accrual_factors.append(tau)
            self.grid_times.append(t)
            prev_dt = next_dt

        #        print(self.grid_times)
        self.accrual_factors = np.array(self.accrual_factors)
        self.num_fwds = len(self.accrual_factors)
        self.fwds = None
        self.use_sobol = None
        self.num_paths = None
        self.numeraire_index = None
        self.fwd_curve = None
        self.vol_curves = None
        self.corr_matrix = None
        self.model_type = None

    #        print("Num FORWARDS", self.num_fwds)

    ###########################################################################

    def simulate_1f(
        self,
        discount_curve,
        vol_curve: IborCapVolCurve,
        num_paths: int = 1000,
        numeraire_index: int = 0,
        use_sobol: bool = True,
        seed: int = 42,
    ):
        """Run the one-factor simulation of the evolution of the forward
        Ibors to generate and store all of the Ibor forward rate paths."""

        if num_paths < 2 or num_paths > 1000000:
            raise FinError("NumPaths must be between 2 and 1 million")

        if discount_curve.value_dt != self.start_dt:
            raise FinError("Curve anchor date not the same as LMM start date.")

        self.num_paths = num_paths
        self.numeraire_index = numeraire_index
        self.use_sobol = use_sobol

        num_grid_points = len(self.grid_dts)

        self.num_fwds = num_grid_points
        self.fwd_curve = []

        for i in range(1, num_grid_points):
            start_dt = self.grid_dts[i - 1]
            end_dt = self.grid_dts[i]
            fwd_rate = discount_curve.fwd_rate(
                start_dt, end_dt, self.float_dc_type
            )
            self.fwd_curve.append(fwd_rate)

        self.fwd_curve = np.array(self.fwd_curve)

        gammas = np.zeros(num_grid_points)
        for ix in range(1, num_grid_points):
            dt = self.grid_dts[ix]
            gammas[ix] = vol_curve.caplet_vol(dt)

        self.fwds = lmm_simulate_fwds_1f(
            self.num_fwds,
            num_paths,
            numeraire_index,
            self.fwd_curve,
            gammas,
            self.accrual_factors,
            use_sobol,
            seed,
        )

    ###########################################################################

    def simulate_mf(
        self,
        discount_curve,
        num_factors: int,
        lambdas: np.ndarray,
        num_paths: int = 10000,
        numeraire_index: int = 0,
        use_sobol: bool = True,
        seed: int = 42,
    ):
        """Run the simulation to generate and store all of the Ibor forward
        rate paths. This is a multi-factorial version so the user must input
        a numpy array consisting of a column for each factor and the number of
        rows must equal the number of grid times on the underlying simulation
        grid. CHECK THIS."""

        #        check_argument_types(self.__init__, locals())

        if num_paths < 2 or num_paths > 1000000:
            raise FinError("NumPaths must be between 2 and 1 million")

        if discount_curve.curve_dt != self.start_dt:
            raise FinError("Curve anchor date not the same as LMM start date.")

        print("LEN LAMBDAS", len(lambdas))
        print("LEN", len(lambdas[0]))
        # We pass a vector of vol discount, one for each factor
        if num_factors != len(lambdas):
            raise FinError("Lambda doesn't have specified number of factors.")

        num_rows = len(lambdas[0])
        if num_rows != self.num_fwds + 1:
            raise FinError("Vol Components needs same number of rows as grid")

        self.num_paths = num_paths
        self.numeraire_index = numeraire_index
        self.use_sobol = use_sobol

        self.num_fwds = len(self.grid_dts) - 1
        self.fwd_curve = []

        for i in range(1, self.num_fwds):
            start_dt = self.grid_dts[i - 1]
            end_dt = self.grid_dts[i]
            fwd_rate = discount_curve.fwd_rate(
                start_dt, end_dt, self.float_dc_type
            )
            self.fwd_curve.append(fwd_rate)

        self.fwd_curve = np.array(self.fwd_curve)

        self.fwds = lmm_simulate_fwds_mf(
            self.num_fwds,
            num_factors,
            num_paths,
            numeraire_index,
            self.fwd_curve,
            lambdas,
            self.accrual_factors,
            use_sobol,
            seed,
        )

    ###########################################################################

    def simulate_nf(
        self,
        discount_curve,
        vol_curve: IborCapVolCurve,
        corr_matrix: np.ndarray,
        model_type: ModelLMMModelTypes,
        num_paths: int = 1000,
        numeraire_index: int = 0,
        use_sobol: bool = True,
        seed: int = 42,
    ):
        """Run the simulation to generate and store all of the Ibor forward
        rate paths using a full factor reduction of the fwd-fwd correlation
        matrix using Cholesky decomposition."""

        check_argument_types(self.__init__, locals())

        if num_paths < 2 or num_paths > 1000000:
            raise FinError("NumPaths must be between 2 and 1 million")

        if isinstance(model_type, ModelLMMModelTypes) is False:
            raise FinError("Model type must be type FinRateModelLMMModelTypes")

        if discount_curve.curve_dt != self.start_dt:
            raise FinError("Curve anchor date not the same as LMM start date.")

        self.num_paths = num_paths
        self.vol_curves = vol_curve
        self.corr_matrix = corr_matrix
        self.model_type = model_type
        self.numeraire_index = numeraire_index
        self.use_sobol = use_sobol

        num_grid_points = len(self.grid_times)

        self.num_fwds = num_grid_points - 1
        self.fwd_curve = []

        for i in range(1, num_grid_points):
            start_dt = self.grid_dts[i - 1]
            end_dt = self.grid_dts[i]
            fwd_rate = discount_curve.forward_rate(
                start_dt, end_dt, self.float_dc_type
            )
            self.fwd_curve.append(fwd_rate)

        self.fwd_curve = np.array(self.fwd_curve)

        zetas = np.zeros(num_grid_points)
        for ix in range(1, num_grid_points):
            dt = self.grid_dts[ix]
            zetas[ix] = vol_curve.caplet_vol(dt)

        # This function does not use Sobol - TODO
        self.fwds = lmm_simulate_fwds_nf(
            self.num_fwds,
            num_paths,
            self.fwd_curve,
            zetas,
            corr_matrix,
            self.accrual_factors,
            seed,
        )

    ###########################################################################

    def value_swaption(
        self,
        settle_dt: Date,
        exercise_dt: Date,
        maturity_dt: Date,
        swaption_type: SwapTypes,
        fixed_cpn: float,
        fixed_freq_type: FrequencyTypes,
        fixed_dc_type: DayCountTypes,
        notional: float = ONE_MILLION,
        float_freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
        float_dc_type: DayCountTypes = DayCountTypes.THIRTY_E_360,
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
    ):
        """Value a swaption in the LMM model using simulated paths of the
        forward curve. This relies on pricing the fixed leg of the swap and
        assuming that the floating leg will be worth par. As a result we only
        need simulate Ibors with the frequency of the fixed leg."""

        # Note that the simulation time steps run all the way out to the last
        # forward rate. However we only really need the forward rates at the
        # expiry date of the option. It may be worth amending the simulate
        # code to impose a limit on the time steps in order to speed up the
        # overall pricing if it requires a new run every time. However once
        # generated, the speed of pricing is not affected so this is not
        # strictly an urgent issue.

        swaption_float_dts = Schedule(
            settle_dt, maturity_dt, float_freq_type, cal_type, bd_type, dg_type
        ).generate()

        for swaption_dt in swaption_float_dts:
            found_dt = False
            for grid_dt in self.grid_dts:
                if swaption_dt == grid_dt:
                    found_dt = True
                    break
            if found_dt is False:
                raise FinError("Swaption float leg not on grid.")

        swaption_fixed_dts = Schedule(
            settle_dt, maturity_dt, fixed_freq_type, cal_type, bd_type, dg_type
        ).generate()

        for swaption_dt in swaption_fixed_dts:
            found_dt = False
            for grid_dt in self.grid_dts:
                if swaption_dt == grid_dt:
                    found_dt = True
                    break
            if found_dt is False:
                raise FinError("Swaption fixed leg not on grid.")

        a = 0
        b = 0

        for grid_dt in self.grid_dts:
            if grid_dt == exercise_dt:
                break
            else:
                a += 1

        for grid_dt in self.grid_dts:
            if grid_dt == maturity_dt:
                break
            else:
                b += 1

        if b == 0:
            raise FinError("Swaption swap maturity date is today.")

        #        num_paths = 1000
        #        v = LMMSwaptionPricer(fixed_cpn, a, b, num_paths,
        #                              fwd0, fwds, taus, is_payer)
        v = 0.0
        return v

    ###########################################################################

    def value_cap_floor(
        self,
        settle_dt: Date,
        maturity_dt: Date,
        cap_floor_type: FinCapFloorTypes,
        cap_floor_rate: float,
        freq_type: FrequencyTypes = FrequencyTypes.QUARTERLY,
        dc_type: DayCountTypes = DayCountTypes.ACT_360,
        notional: float = ONE_MILLION,
        cal_type: CalendarTypes = CalendarTypes.WEEKEND,
        bd_type: BusDayAdjustTypes = BusDayAdjustTypes.FOLLOWING,
        dg_type: DateGenRuleTypes = DateGenRuleTypes.BACKWARD,
    ):
        """Value a cap or floor in the LMM."""

        cap_floor_dts = Schedule(
            settle_dt, maturity_dt, freq_type, cal_type, bd_type, dg_type
        ).generate()

        for cap_floorlet_dt in cap_floor_dts:
            found_dt = False
            for grid_dt in self.grid_dts:
                if cap_floorlet_dt == grid_dt:
                    found_dt = True
                    break
            if found_dt is False:
                raise FinError("CapFloor date not on grid.")

        num_fwds = len(cap_floor_dts)
        num_paths = self.num_paths
        K = cap_floor_rate

        is_cap = 0
        if cap_floor_type == FinCapFloorTypes.CAP:
            is_cap = 1

        fwd0 = self.fwd_curve
        fwds = self.fwds
        taus = self.accrual_factors

        v = lmm_cap_flr_pricer(
            num_fwds, num_paths, K, fwd0, fwds, taus, is_cap
        )

        # Sum the cap/floorlets to get cap/floor value
        v_cap_floor = 0.0
        for v_cap_floor_let in v:
            v_cap_floor += v_cap_floor_let * notional

        return v_cap_floor

    ###########################################################################

    def __repr__(self):
        """Function to allow us to print the LMM Products details."""

        s = "Function not written"
        return s

    ###########################################################################

    def _print(self):
        """Alternative print method."""

        print(self)


###############################################################################
