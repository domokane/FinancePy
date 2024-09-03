import copy
import numpy as np
import pandas as pd
from scipy import optimize


from ...utils.date import datediff
from ...utils.global_vars import g_days_in_year
from ...products.rates.ibor_single_curve import IborSingleCurve
from ...products.rates.ibor_benchmarks_report import ibor_benchmarks_report


class IborSingleCurveSmoothingCalibrator(object):
    """
    Non-parametric fitting of a curve with smoothness. We use dfs
    for all coupon (cashflow) dates as input variables and impose a smoothness
    penalty on (at the moment) an approximation to the second derivative of the
    yields (zero rates). With non-zero smoothness penalty the fit to
    the benchmarks is not exact.
    """

    def __init__(self, ibor_curve: IborSingleCurve):
        """
        Initialize with an (unbuilt) ibor curve. Note that this curve is not modified
        during fitting as we take a deep copy
        """

        #
        self._curve = copy.deepcopy(ibor_curve)

        self._collect_all_knot_dates()

    def _collect_all_knot_dates(self):
        """
        Collect all dats on which the discount factors are explicitly required
        to price bechmarks. Interpolation only used for times in between knot dates
        """

        # shorter name
        c = self._curve

        dates = []
        dates.append(self._curve.value_dt)

        for depo in c.used_deposits:
            dates.append(depo.start_dt)
            dates.append(depo.maturity_dt)

        for fra in c.used_fras:
            dates.append(fra.start_dt)
            dates.append(fra.maturity_dt)

        for swap in c.used_swaps:
            dates.extend(swap.fixed_leg.payment_dts)
            dates.append(swap.maturity_dt)

        # remove duplicates
        dates = list(set(dates))
        dates.sort()
        self._knot_dates = dates
        self._knot_times = np.array(
            [(d - dates[0]) / g_days_in_year for d in dates]
        )

    def _repricing_objectives(self, curve_to_use=None):

        if curve_to_use is not None:
            curve = curve_to_use
        else:
            curve = self._curve

        valuation_date = curve.value_dt
        out = np.zeros(
            len(curve.used_deposits)
            + len(curve.used_fras)
            + len(curve.used_swaps)
        )

        idx = 0
        for depo in curve.used_deposits:
            # do not need to be too exact here
            acc_factor = datediff(depo.start_dt, depo.maturity_dt)
            # as rate
            r = (
                -np.log(depo.value(valuation_date, curve) / depo.notional)
                / acc_factor
            )
            out[idx] = r
            idx = idx + 1

        for fra in curve.used_fras:
            # do not need to be too exact here
            acc_factor = datediff(fra.start_dt, fra.maturity_dt)
            v = fra.value(valuation_date, curve) / fra.notional / acc_factor
            out[idx] = v
            idx = idx + 1

        for swap in curve.used_swaps:
            v = (
                swap.value(valuation_date, curve)
                / swap.fixed_leg.notional
                / swap.pv01(valuation_date, curve)
            )
            out[idx] = v
            idx = idx + 1

        return out

    def _smoothing_objectives(self, curve_to_use=None):
        if curve_to_use is not None:
            curve = curve_to_use
        else:
            curve = self._curve

        start_times = curve._times[:-1]
        end_times = curve._times[1:]
        tenors = end_times - start_times

        start_dfs = curve._dfs[:-1]
        end_dfs = curve._dfs[1:]
        fdfs = end_dfs / start_dfs

        # forward cc rates -- first derivative of the yields
        fcc_rates = -np.log(fdfs) / tenors
        fcc_rate_derivs = (
            np.diff(fcc_rates) / tenors[1:]
        )  # some qs about which tenor we should divide this by

        return fcc_rate_derivs

    def fit(self, smoothness=1e-6, report_progress=False) -> IborSingleCurve:
        """
        fit the curve with a given smoothness
        """

        def _obj_f(dfs):
            curve = self._curve
            curve._dfs[1:] = dfs

            curve._interpolator.fit(curve._times, curve._dfs)

            fit_tgts = self._repricing_objectives(curve)
            smth_tgts = self._smoothing_objectives(curve)
            smth_tgts *= smoothness

            if report_progress:
                n_fit = len(fit_tgts)
                n_smth = len(smth_tgts)
                print(
                    f"fit rmse in bps={10000*np.linalg.norm(fit_tgts)/n_fit}, smooth = {10000*np.linalg.norm(smth_tgts)/n_smth}"
                )

            return np.concatenate((fit_tgts, smth_tgts))

        # If the curve passed during construction is already built, use that
        # for the initial guess of our optimization. If not, use the bootstrap
        if self._curve._is_built:
            init_curve = self._curve
        else:
            init_curve = copy.deepcopy(self._curve)
            init_curve._check_refit = False
            init_curve._build_curve_using_1d_solver(
                **init_curve._optional_interp_params
            )

        dfs0 = init_curve.df(self._knot_dates)

        self._curve._times = self._knot_times
        self._curve._dfs = np.ones_like(self._knot_times, dtype=float)
        self._curve._is_built = True

        res = optimize.least_squares(
            _obj_f, dfs0[1:], bounds=(0, np.inf), ftol=1e-4, xtol=1e-6
        )

        self._curve._dfs[1:] = np.array(res.x)
        self._curve._interpolator.fit(self._curve._times, self._curve._dfs)

        fit_report = self._generate_fit_report(smoothness)

        return copy.deepcopy(self._curve), fit_report

    def _generate_fit_report(self, smoothness, curve_to_use=None):
        if curve_to_use is not None:
            curve = curve_to_use
        else:
            curve = self._curve

        fit_tgts = self._repricing_objectives(curve)
        smth_tgts = self._smoothing_objectives(curve)

        fit_labels = [f"fit_{n:03d}" for n in np.arange(len(fit_tgts))]
        smth_labels = [f"smth_{n:03d}" for n in np.arange(len(smth_tgts))]

        df1 = pd.DataFrame(
            columns=["tgt_label", "value_in_bps"],
            data=zip(fit_labels, 10000 * fit_tgts),
        )
        df2 = pd.DataFrame(
            columns=["tgt_label", "value_in_bps"],
            data=zip(smth_labels, 10000 * smth_tgts),
        )

        # decorate df1 with extra info
        df_bmi = ibor_benchmarks_report(curve)
        df_bmi["tgt_label"] = df1["tgt_label"]
        df1 = df1.merge(df_bmi, how="inner", on="tgt_label")

        # decorate df2 with extra info
        df2["type"] = "d2yield_dt2"
        df2["start_date"] = self._knot_dates[1:-1]
        df2["maturity_dt"] = self._knot_dates[2:]

        df = pd.concat((df1, df2), axis=0, sort=False).reset_index(drop=True)
        df["smoothness"] = smoothness
        return df
