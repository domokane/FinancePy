# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


from .model import Model
from ..utils.helpers import check_argument_types
from typing import Union
import numpy as np
from ..utils.global_types import OptionTypes
from .bs_american_analytic import BlackScholesAlgorithmTypes
from .bs_american_analytic import american_value
from .bs_european_analytic import european_value
from ..utils.error import FinError

# This file defines a flexible Black-Scholes pricing engine class from the
# user's FinancePy library. It supports multiple pricing methodologies for
# both European and American equity options.

# UNDER CONSTRUCTION - REDESIGN IN PROGRESS!!!!!!!!!!!!

########################################################################################


class BlackScholes(Model):
    """
    Black-Scholes model class supporting several pricing methods.
    """

    def __init__(
        self,
        volatility: Union[float, np.ndarray],
        # bs_type: BlackScholesAlgorithmTypes = BlackScholesAlgorithmTypes.DEFAULT,
        # num_steps_per_year: int = 52,
        # num_paths: int = 10000,
        # seed: int = 42,
        # use_sobol: bool = False,
        # params: dict | None = None,
    ):
        check_argument_types(self.__init__, locals())

        self.volatility = volatility
#        self.bs_type = bs_type
#        self.num_steps_per_year = num_steps_per_year
#        self.num_paths = num_paths
#        self.seed = seed
#        self.use_sobol = use_sobol
#        self.params = params or {}

#        self.poly_degree = self.params.get("poly_degree", 3)
#        self.fit_type = self.params.get("fit_type", BoundaryFitTypes.HERMITE_E)

    ####################################################################################

    def value(
        self,
        spot_price: float,
        time_to_expiry: float,
        strike_price: float,
        risk_free_rate: float,
        dividend_rate: float,
        opt_type: OptionTypes,
        bs_type: BlackScholesAlgorithmTypes
    ) -> float:
        """
        Compute the option value using the selected Black-Scholes method.
        """

        if opt_type in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
            bs_type = (
                BlackScholesAlgorithmTypes.ANALYTICAL
                if bs_type is BlackScholesAlgorithmTypes.DEFAULT
                else bs_type
            )
            return european_value(
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
                BlackScholesAlgorithmTypes.CRR_TREE
                if bs_type is BlackScholesAlgorithmTypes.DEFAULT
                else bs_type
            )
            return american_value(
                bs_type,
                spot_price,
                time_to_expiry,
                strike_price,
                risk_free_rate,
                dividend_rate,
                opt_type,
            )

        raise FinError(f"Unsupported option type: {opt_type}")

    # ####################################################################################

    # def _value_european(
    #     self,
    #     bs_type: BlackScholesAlgorithmTypes,
    #     spot_price: float,
    #     time_to_expiry: float,
    #     strike_price: float,
    #     risk_free_rate: float,
    #     dividend_rate: float,
    #     opt_type: OptionTypes,
    # ) -> float:

    #     if bs_type is BlackScholesAlgorithmTypes.ANALYTICAL:
    #         return european_value(
    #             spot_price,
    #             time_to_expiry,
    #             strike_price,
    #             risk_free_rate,
    #             dividend_rate,
    #             self.volatility,
    #             opt_type.value,
    #         )

    #     if bs_type is BlackScholesAlgorithmTypes.CRR_TREE:
    #         return self._value_crr_tree(
    #             spot_price,
    #             time_to_expiry,
    #             strike_price,
    #             risk_free_rate,
    #             dividend_rate,
    #             opt_type,
    #         )

    #     if bs_type is BlackScholesAlgorithmTypes.LSMC:
    #         return self._value_lsmc(
    #             spot_price,
    #             time_to_expiry,
    #             strike_price,
    #             risk_free_rate,
    #             dividend_rate,
    #             opt_type,
    #         )

    #     if bs_type is BlackScholesAlgorithmTypes.FINITE_DIFFERENCE:
    #         return self._value_fd(
    #             spot_price,
    #             time_to_expiry,
    #             strike_price,
    #             risk_free_rate,
    #             dividend_rate,
    #             opt_type,
    #         )

    #     if bs_type is BlackScholesAlgorithmTypes.PSOR:
    #         return self._value_psor(
    #             spot_price,
    #             time_to_expiry,
    #             strike_price,
    #             risk_free_rate,
    #             dividend_rate,
    #             opt_type,
    #         )

    #     raise FinError(f"Unsupported European Black-Scholes type: {bs_type}")

    # ####################################################################################

    # def _value_american(
    #     self,
    #     bs_type: BlackScholesAlgorithmTypes,
    #     spot_price: float,
    #     time_to_expiry: float,
    #     strike_price: float,
    #     risk_free_rate: float,
    #     dividend_rate: float,
    #     opt_type: OptionTypes,
    # ) -> float:

    #     if bs_type is BlackScholesAlgorithmTypes.BARONE_ADESI:
    #         return baw_value(
    #             spot_price,
    #             time_to_expiry,
    #             strike_price,
    #             risk_free_rate,
    #             dividend_rate,
    #             self.volatility,
    #             opt_type.value,
    #         )

    #     if bs_type is BlackScholesAlgorithmTypes.BJERKSUND_STENSLAND:
    #         return bjerksund_stensland_value(
    #             spot_price,
    #             time_to_expiry,
    #             strike_price,
    #             risk_free_rate,
    #             dividend_rate,
    #             self.volatility,
    #             opt_type.value,
    #         )

    #     if bs_type is BlackScholesAlgorithmTypes.CRR_TREE:
    #         return self._value_crr_tree(
    #             spot_price,
    #             time_to_expiry,
    #             strike_price,
    #             risk_free_rate,
    #             dividend_rate,
    #             opt_type,
    #         )

    #     if bs_type is BlackScholesAlgorithmTypes.LSMC:
    #         return self._value_lsmc(
    #             spot_price,
    #             time_to_expiry,
    #             strike_price,
    #             risk_free_rate,
    #             dividend_rate,
    #             opt_type,
    #         )

    #     if bs_type is BlackScholesAlgorithmTypes.FINITE_DIFFERENCE:
    #         return self._value_fd(
    #             spot_price,
    #             time_to_expiry,
    #             strike_price,
    #             risk_free_rate,
    #             dividend_rate,
    #             opt_type,
    #         )

    #     if bs_type is BlackScholesAlgorithmTypes.PSOR:
    #         return self._value_psor(
    #             spot_price,
    #             time_to_expiry,
    #             strike_price,
    #             risk_free_rate,
    #             dividend_rate,
    #             opt_type,
    #         )

    #     raise FinError(f"Unsupported American Black-Scholes type: {bs_type}")

    # ####################################################################################

    # def _value_crr_tree(
    #     self,
    #     spot_price,
    #     time_to_expiry,
    #     strike_price,
    #     risk_free_rate,
    #     dividend_rate,
    #     opt_type,
    # ) -> float:
    #     return crr_tree_val_avg(
    #         spot_price,
    #         risk_free_rate,
    #         dividend_rate,
    #         self.volatility,
    #         self.num_steps_per_year,
    #         time_to_expiry,
    #         opt_type.value,
    #         strike_price,
    #     )["value"]

    # ####################################################################################

    # def _value_lsmc(
    #     self,
    #     spot_price,
    #     time_to_expiry,
    #     strike_price,
    #     risk_free_rate,
    #     dividend_rate,
    #     opt_type,
    # ) -> float:
    #     return equity_lsmc(
    #         spot_price=spot_price,
    #         risk_free_rate=risk_free_rate,
    #         dividend_yield=dividend_rate,
    #         sigma=self.volatility,
    #         num_steps_per_year=self.num_steps_per_year,
    #         num_paths=self.num_paths,
    #         time_to_expiry=time_to_expiry,
    #         opt_type_value=opt_type.value,
    #         strike_price=strike_price,
    #         poly_degree=self.poly_degree,
    #         fit_type_value=self.fit_type.value,
    #         use_sobol=self.use_sobol,
    #         seed=self.seed,
    #     )

    # ####################################################################################

    # def _value_fd(
    #     self,
    #     spot_price,
    #     time_to_expiry,
    #     strike_price,
    #     risk_free_rate,
    #     dividend_rate,
    #     opt_type,
    # ) -> float:
    #     return black_scholes_fd(
    #         spot_price=spot_price,
    #         time_to_expiry=time_to_expiry,
    #         strike_price=strike_price,
    #         risk_free_rate=risk_free_rate,
    #         dividend_yield=dividend_rate,
    #         volatility=self.volatility,
    #         opt_type=opt_type.value,
    #         **self.params,
    #     )

    # ####################################################################################

    # def _value_psor(
    #     self,
    #     spot_price,
    #     time_to_expiry,
    #     strike_price,
    #     risk_free_rate,
    #     dividend_rate,
    #     opt_type,
    # ) -> float:
    #     return black_scholes_fd_psor(
    #         spot_price=spot_price,
    #         time_to_expiry=time_to_expiry,
    #         strike_price=strike_price,
    #         risk_free_rate=risk_free_rate,
    #         dividend_yield=dividend_rate,
    #         volatility=self.volatility,
    #         opt_type=opt_type.value,
    #         **self.params,
    #     )
