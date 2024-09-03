##############################################################################
# Copyright (C) 2018, 2019, 2020 Dominic O'Kane
##############################################################################


from ...models.black_scholes import BlackScholes
from ...utils.global_vars import g_days_in_year
from ...utils.date import Date

##########################################################################

bump = 1e-4

##########################################################################


class FXOption:
    """Class that is used to perform perturbation risk for FX options."""

    def value(
        self,
        value_dt: Date,
        spot_fx_rate: float,
        domestic_curve,
        foreign_curve,
        model,
    ):

        print("You should not be here!")
        return 0.0

    ###########################################################################

    def delta(
        self, value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    ):
        """Calculate the option delta (FX rate sensitivity) by adding on a
        small bump and calculating the change in the option price."""

        v = self.value(
            value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
        )

        v_bumped = self.value(
            value_dt, spot_fx_rate + bump, domestic_curve, foreign_curve, model
        )

        if isinstance(v_bumped, dict):
            delta = (v_bumped["value"] - v["value"]) / bump
        else:
            delta = (v_bumped - v) / bump

        return delta

    ###########################################################################

    def gamma(
        self, value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    ):
        """Calculate the option gamma (delta sensitivity) by adding on a
        small bump and calculating the change in the option delta."""

        v = self.delta(
            value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
        )

        v_bumped_dn = self.delta(
            value_dt, spot_fx_rate + bump, domestic_curve, foreign_curve, model
        )

        v_bumped_up = self.delta(
            value_dt, spot_fx_rate + bump, domestic_curve, foreign_curve, model
        )

        if isinstance(v, dict):
            num = (
                v_bumped_up["value"] - 2.0 * v["value"] + v_bumped_dn["value"]
            )
            gamma = num / bump / 2.0
        else:
            gamma = (v_bumped_up - 2.0 * v + v_bumped_dn) / bump / 2.0

        return gamma

    ###########################################################################

    def vega(
        self, value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    ):
        """Calculate the option vega (volatility sensitivity) by adding on a
        small bump and calculating the change in the option price."""

        bump = 0.01

        v = self.value(
            value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
        )

        vp = self.value(
            value_dt,
            spot_fx_rate,
            domestic_curve,
            foreign_curve,
            BlackScholes(model.volatility + bump),
        )

        if isinstance(v, dict):
            vega = vp["value"] - v["value"]  # / bump
        else:
            vega = vp - v  # / bump

        return vega

    ###########################################################################

    def theta(
        self, value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    ):
        """Calculate the option theta (calendar time sensitivity) by moving
        forward one day and calculating the change in the option price."""

        v = self.value(
            value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
        )

        next_dt = value_dt.add_days(1)

        domestic_curve.value_dt = next_dt
        foreign_curve.value_dt = next_dt

        v_bumped = self.value(
            next_dt, spot_fx_rate, domestic_curve, foreign_curve, model
        )

        bump = 1.0 / g_days_in_year

        if isinstance(v, dict):
            theta = (v_bumped["value"] - v["value"]) / bump
        else:
            theta = (v_bumped - v) / bump

        # Don't forget to reset the value dates
        domestic_curve.value_dt = value_dt
        foreign_curve.value_dt = value_dt

        return theta

    ###########################################################################

    def rho(
        self, value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
    ):
        """Calculate the option rho (interest rate sensitivity) by perturbing
        the discount curve and revaluing."""

        v = self.value(
            value_dt, spot_fx_rate, domestic_curve, foreign_curve, model
        )
        v_bumped = self.value(
            value_dt,
            spot_fx_rate,
            domestic_curve.bump(bump),
            foreign_curve,
            model,
        )

        if isinstance(v, dict):
            rho = (v_bumped["value"] - v["value"]) / bump
        else:
            rho = (v_bumped - v) / bump

        return rho


##########################################################################
##########################################################################
