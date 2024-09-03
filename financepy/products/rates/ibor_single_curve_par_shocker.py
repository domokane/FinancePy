import copy
from typing import Union
import numpy as np

from ...utils.global_vars import g_basis_point
from ...products.rates.ibor_deposit import IborDeposit
from ...products.rates.ibor_fra import IborFRA
from ...products.rates.ibor_swap import IborSwap
from ...products.rates.ibor_single_curve import IborSingleCurve
from ...products.rates.ibor_benchmarks_report import (
    benchmarks_report,
    ibor_benchmarks_report,
)


class IborSingleCurveParShocker:
    """
    A class to apply par-rate, ie benchmark, bumps to a Libor curve. Takes a base curve
    and provides methods to apply bumps that return bumped curves
    """

    def __init__(self, base_curve: IborSingleCurve):
        """
        Init with a base curve to be bumped. Bumps do not affect this curve.

        """
        self._base_curve = base_curve
        self._benchmarks_report = ibor_benchmarks_report(
            self._base_curve, include_objects=True
        )

    def benchmarks_report(self):
        """
        Access the benchmarks report that we create when the shocker is initialized
        """
        return self._benchmarks_report

    def n_benchmarks(self):
        """
        Total number of benchmarks
        """
        return len(self._benchmarks_report)

    def apply_bump_to_benchmark(
        self, benchmark_idx: int, bump_size=1.0 * g_basis_point
    ):
        """
        Apply a shock of a given size to a given bechmark. Indexing is per the benchmark report
        """
        composite_shock = np.zeros(self.n_benchmarks())
        composite_shock[benchmark_idx] = bump_size
        return self.apply_composite_bump(composite_shock)

    def apply_composite_bump(self, bump_sizes: Union[np.array, list]):
        """Apply a composite bump to base_curve. A composite bump is a list/array of bumps, one per bechmark

        Args:
            bump_sizes (Union[np.array, list]): a list/array of bump sizes, one per benchmark

        Returns:
            IborSingleCurve: A bumped curve
        """
        bumped_depos = []
        bumped_fras = []
        bumped_swaps = []

        for benchmark, bump_size in zip(
            self._benchmarks_report["benchmark_objects"].values, bump_sizes
        ):
            bumped_benchmark = copy.deepcopy(benchmark)
            if isinstance(bumped_benchmark, IborDeposit):
                bumped_benchmark.deposit_rate += bump_size
                bumped_depos.append(bumped_benchmark)
            if isinstance(bumped_benchmark, IborFRA):
                bumped_benchmark.fra_rate += bump_size
                bumped_fras.append(bumped_benchmark)
            if isinstance(bumped_benchmark, IborSwap):
                bumped_benchmark.set_fixed_rate(
                    bumped_benchmark.get_fixed_rate() + bump_size
                )
                bumped_swaps.append(bumped_benchmark)

        # This assumes that the base curve was built using the default method as defined
        # in IborSingleCurve._build_curve() based on interp_type. This at the moment excludes
        # the non-parametric smoothing calibration (for which, incidentally, par bumps are not a giid
        # idea anyway). But in the future we should keep track of what exactly we used to build
        # the base curve and use the same method here
        bumped_curve = IborSingleCurve(
            self._base_curve.value_dt,
            bumped_depos,
            bumped_fras,
            bumped_swaps,
            interp_type=self._base_curve._interp_type,
            check_refit=self._base_curve._check_refit,
            do_build=True,
            **self._base_curve._optional_interp_params
        )

        return bumped_curve
