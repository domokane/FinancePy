# Copyright (C) 2018, 2019, 2020 Dominic O'Kane

from typing import Union
from enum import Enum
import numpy as np

from ..utils.global_types import OptionTypes
from ..utils.error import FinError
from ..utils.helpers import check_argument_types
from .model import Model
from .equity_crr_tree import crr_tree_val_avg
from .equity_lsmc import equity_lsmc, BoundaryFitTypes

# This file defines a flexible Black-Scholes pricing engine class from the
# user's FinancePy library. It supports multiple pricing methodologies for
# both European and American equity options.

from .black_scholes_analytic import (
    bs_value,
    baw_value,
    bjerksund_stensland_value,
)


from .finite_difference import black_scholes_fd
from .finite_difference_psor import black_scholes_fd_psor

########################################################################################


class BlackScholesTypes(Enum):

    DEFAULT = 0
    ANALYTICAL = 1
    CRR_TREE = 2
    BARONE_ADESI = 3
    LSMC = 4
    BJERKSUND_STENSLAND = 5
    FINITE_DIFFERENCE = 6
    PSOR = 7


########################################################################################


class BlackScholes(Model):
    """
    Black-Scholes model class supporting several pricing methods.
    """

    def __init__(
        self,
        volatility: Union[float, np.ndarray],
        bs_type: BlackScholesTypes = BlackScholesTypes.DEFAULT,
        num_steps_per_year: int = 52,
        num_paths: int = 10000,
        seed: int = 42,
        use_sobol: bool = False,
        params: dict | None = None,
    ):
        check_argument_types(self.__init__, locals())

        self.volatility = volatility
        self.bs_type = bs_type
        self.num_steps_per_year = num_steps_per_year
        self.num_paths = num_paths
        self.seed = seed
        self.use_sobol = use_sobol
        self.params = params or {}

        self.poly_degree = self.params.get("poly_degree", 3)
        self.fit_type = self.params.get("fit_type", BoundaryFitTypes.HERMITE_E)

    def value(
        self,
        spot_price: float,
        time_to_expiry: float,
        strike_price: float,
        risk_free_rate: float,
        dividend_rate: float,
        opt_type: OptionTypes,
    ) -> float:
        """
        Compute the option value using the selected Black-Scholes method.
        """

        if opt_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
            bs_type = (
                BlackScholesTypes.ANALYTICAL
                if self.bs_type is BlackScholesTypes.DEFAULT
                else self.bs_type
            )
            return self._value_european(
                bs_type,
                spot_price,
                time_to_expiry,
                strike_price,
                risk_free_rate,
                dividend_rate,
                opt_type,
            )

        if opt_type in (OptionTypes.AMERICAN_CALL, OptionTypes.AMERICAN_PUT):
            bs_type = (
                BlackScholesTypes.CRR_TREE
                if self.bs_type is BlackScholesTypes.DEFAULT
                else self.bs_type
            )
            return self._value_american(
                bs_type,
                spot_price,
                time_to_expiry,
                strike_price,
                risk_free_rate,
                dividend_rate,
                opt_type,
            )

        raise FinError(f"Unsupported option type: {opt_type}")

    def _value_european(
        self,
        bs_type: BlackScholesTypes,
        spot_price: float,
        time_to_expiry: float,
        strike_price: float,
        risk_free_rate: float,
        dividend_rate: float,
        opt_type: OptionTypes,
    ) -> float:

        if bs_type is BlackScholesTypes.ANALYTICAL:
            return bs_value(
                spot_price,
                time_to_expiry,
                strike_price,
                risk_free_rate,
                dividend_rate,
                self.volatility,
                opt_type.value,
            )

        if bs_type is BlackScholesTypes.CRR_TREE:
            return self._value_crr_tree(
                spot_price,
                time_to_expiry,
                strike_price,
                risk_free_rate,
                dividend_rate,
                opt_type,
            )

        if bs_type is BlackScholesTypes.LSMC:
            return self._value_lsmc(
                spot_price,
                time_to_expiry,
                strike_price,
                risk_free_rate,
                dividend_rate,
                opt_type,
            )

        if bs_type is BlackScholesTypes.FINITE_DIFFERENCE:
            return self._value_fd(
                spot_price,
                time_to_expiry,
                strike_price,
                risk_free_rate,
                dividend_rate,
                opt_type,
            )

        if bs_type is BlackScholesTypes.PSOR:
            return self._value_psor(
                spot_price,
                time_to_expiry,
                strike_price,
                risk_free_rate,
                dividend_rate,
                opt_type,
            )

        raise FinError(f"Unsupported European Black-Scholes type: {bs_type}")

    def _value_american(
        self,
        bs_type: BlackScholesTypes,
        spot_price: float,
        time_to_expiry: float,
        strike_price: float,
        risk_free_rate: float,
        dividend_rate: float,
        opt_type: OptionTypes,
    ) -> float:

        if bs_type is BlackScholesTypes.BARONE_ADESI:
            phi = 1 if opt_type is OptionTypes.AMERICAN_CALL else -1
            return baw_value(
                spot_price,
                time_to_expiry,
                strike_price,
                risk_free_rate,
                dividend_rate,
                self.volatility,
                phi,
            )

        if bs_type is BlackScholesTypes.BJERKSUND_STENSLAND:
            return bjerksund_stensland_value(
                spot_price,
                time_to_expiry,
                strike_price,
                risk_free_rate,
                dividend_rate,
                self.volatility,
                opt_type.value,
            )

        if bs_type is BlackScholesTypes.CRR_TREE:
            return self._value_crr_tree(
                spot_price,
                time_to_expiry,
                strike_price,
                risk_free_rate,
                dividend_rate,
                opt_type,
            )

        if bs_type is BlackScholesTypes.LSMC:
            return self._value_lsmc(
                spot_price,
                time_to_expiry,
                strike_price,
                risk_free_rate,
                dividend_rate,
                opt_type,
            )

        if bs_type is BlackScholesTypes.FINITE_DIFFERENCE:
            return self._value_fd(
                spot_price,
                time_to_expiry,
                strike_price,
                risk_free_rate,
                dividend_rate,
                opt_type,
            )

        if bs_type is BlackScholesTypes.PSOR:
            return self._value_psor(
                spot_price,
                time_to_expiry,
                strike_price,
                risk_free_rate,
                dividend_rate,
                opt_type,
            )

        raise FinError(f"Unsupported American Black-Scholes type: {bs_type}")

    def _value_crr_tree(
        self,
        spot_price,
        time_to_expiry,
        strike_price,
        risk_free_rate,
        dividend_rate,
        opt_type,
    ) -> float:
        return crr_tree_val_avg(
            spot_price,
            risk_free_rate,
            dividend_rate,
            self.volatility,
            self.num_steps_per_year,
            time_to_expiry,
            opt_type.value,
            strike_price,
        )["value"]

    def _value_lsmc(
        self,
        spot_price,
        time_to_expiry,
        strike_price,
        risk_free_rate,
        dividend_rate,
        opt_type,
    ) -> float:
        return equity_lsmc(
            spot_price=spot_price,
            risk_free_rate=risk_free_rate,
            dividend_yield=dividend_rate,
            sigma=self.volatility,
            num_steps_per_year=self.num_steps_per_year,
            num_paths=self.num_paths,
            time_to_expiry=time_to_expiry,
            opt_type_value=opt_type.value,
            strike_price=strike_price,
            poly_degree=self.poly_degree,
            fit_type_value=self.fit_type.value,
            use_sobol=self.use_sobol,
            seed=self.seed,
        )

    def _value_fd(
        self,
        spot_price,
        time_to_expiry,
        strike_price,
        risk_free_rate,
        dividend_rate,
        opt_type,
    ) -> float:
        return black_scholes_fd(
            spot_price=spot_price,
            time_to_expiry=time_to_expiry,
            strike_price=strike_price,
            risk_free_rate=risk_free_rate,
            dividend_yield=dividend_rate,
            volatility=self.volatility,
            opt_type=opt_type.value,
            **self.params,
        )

    def _value_psor(
        self,
        spot_price,
        time_to_expiry,
        strike_price,
        risk_free_rate,
        dividend_rate,
        opt_type,
    ) -> float:
        return black_scholes_fd_psor(
            spot_price=spot_price,
            time_to_expiry=time_to_expiry,
            strike_price=strike_price,
            risk_free_rate=risk_free_rate,
            dividend_yield=dividend_rate,
            volatility=self.volatility,
            opt_type=opt_type.value,
            **self.params,
        )
