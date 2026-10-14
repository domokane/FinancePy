##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from ...utils.error import FinError
from ...utils.global_types import OptionTypes

from ...utils.helpers import label_to_string, check_argument_types
from ...utils.date import Date
from ...products.equity.equity_option import EquityOption
from ...market.curves.flat_discount_curve import DiscountCurve
from ...utils.check_values import check_curve_dt
from ...utils.check_values import check_stock_price
from ...utils.helpers import option_years

from ...models.bs_european_analytic import european_value
from ...models.black_scholes import BlackScholes
from ...models.model import Model

########################################################################################
# TODO: Do we need to day count adjust option payoffs ?
# TODO: Monte Carlo pricer
########################################################################################


class EquityForwardStartOption(EquityOption):
    """European forward-start equity option.

    The strike is set equal to the stock price on start_dt.
    The option expires on final_expiry_dt.

    The option is valued under Black-Scholes using forward interest
    and dividend rates between the strike-setting and expiry dates.
    """

    def __init__(
        self,
        start_dt: Date,
        final_expiry_dt: Date,
        opt_type: OptionTypes,
    ):
        """Create the EquityForwardStartOption by passing in the start date on
        which the first strike is set, the final expiry date and whether it is a
        call or a put. The frequency determines the reset dates between the two
        dates; the calendar and business day conventions adjust those dates."""

        check_argument_types(self.__init__, locals())

        if opt_type not in (OptionTypes.EUROPEAN_CALL, OptionTypes.EUROPEAN_PUT):
            raise FinError("Unknown option type")

        if final_expiry_dt <= start_dt:
            raise FinError("Expiry date must be after start date")

        self.start_dt = start_dt
        self.final_expiry_dt = final_expiry_dt
        self.opt_type = opt_type

    ###########################################################################

    def value(
        self,
        value_dt: Date,
        stock_price: float,
        discount_curve: DiscountCurve,
        dividend_curve: DiscountCurve,
        model: Model,
    ):
        if not isinstance(value_dt, Date):
            raise FinError("Valuation date is not a Date")

        check_curve_dt(value_dt, discount_curve)
        check_curve_dt(value_dt, dividend_curve)

        if value_dt > self.start_dt:
            raise FinError(
                "Value date after strike-setting date"
            )

        if not isinstance(model, BlackScholes):
            raise FinError("Unknown model type")

        check_stock_price(stock_price)

        t = option_years(self.start_dt, self.final_expiry_dt)

        fwd_r = discount_curve.fwd_zero_rate_cc(
            self.start_dt, self.final_expiry_dt
        )

        fwd_q = dividend_curve.fwd_zero_rate_cc(
            self.start_dt, self.final_expiry_dt
        )

        dq_start = dividend_curve.df(self.start_dt)

        volatility = max(model.volatility, 1e-6)

        v = european_value(
            1.0,
            t,
            1.0,
            fwd_r,
            fwd_q,
            volatility,
            self.opt_type.value,
        )

        return stock_price * dq_start * v

    ###########################################################################

    def __repr__(self):
        s = label_to_string("OBJECT_TYPE", type(self).__name__)
        s += label_to_string("START_DATE", self.start_dt)
        s += label_to_string("FINAL EXPIRY DATE", self.final_expiry_dt)
        s += label_to_string("OPTION_TYPE", self.opt_type)
        return s

    ###########################################################################

    def _print(self):
        """Simple print function for backward compatibility."""
        print(self)
