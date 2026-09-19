# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 10:58:29 2026

@author: Dominic
"""

# volatility_surface.py

import numpy as np

from ..utils.error import FinError

########################################################################################


class ImpliedVolatilitySurface:
    """
    Base class for implied-volatility surface models.

    Concrete models should implement

        implied_volatility(forward, strike, t_exp)

    and may optionally implement

        total_variance(forward, strike, t_exp).

    Examples include

        SVI
        SSVI
        SABR
        interpolated market surfaces
        stochastic-volatility implied-volatility surfaces
    """

    def implied_volatility(
        self,
        forward,
        strike,
        t_exp,
    ):
        raise NotImplementedError

    ####################################################################################

    def total_variance(
        self,
        forward,
        strike,
        t_exp,
    ):
        vol = self.implied_volatility(
            forward,
            strike,
            t_exp,
        )

        return vol * vol * t_exp

    ####################################################################################

    def implied_volatility_curve(
        self,
        forward,
        strikes,
        t_exp,
    ):

        strikes = np.asarray(
            strikes,
            dtype=float,
        )

        if strikes.ndim != 1:
            raise FinError("Strikes must be one-dimensional.")

        if forward <= 0.0:
            raise FinError("Forward must be positive.")

        if t_exp <= 0.0:
            raise FinError("Expiry must be positive.")

        vols = np.empty(
            len(strikes),
            dtype=float,
        )

        for i, strike in enumerate(strikes):

            vols[i] = self.implied_volatility(
                forward,
                strike,
                t_exp,
            )

        return vols

    ####################################################################################

    def implied_volatility_surface(
        self,
        forwards,
        strikes,
        expiries,
    ):

        forwards = np.asarray(
            forwards,
            dtype=float,
        )

        strikes = np.asarray(
            strikes,
            dtype=float,
        )

        expiries = np.asarray(
            expiries,
            dtype=float,
        )

        if expiries.ndim != 1:
            raise FinError("Expiries must be one-dimensional.")

        if strikes.ndim != 1:
            raise FinError("Strikes must be one-dimensional.")

        if forwards.ndim != 1:
            raise FinError("Forwards must be one-dimensional.")

        if len(forwards) != len(expiries):
            raise FinError("Number of forwards must equal number of expiries.")

        vols = np.empty(
            (
                len(expiries),
                len(strikes),
            ),
            dtype=float,
        )

        for i, t_exp in enumerate(expiries):

            for j, strike in enumerate(strikes):

                vols[i, j] = self.implied_volatility(
                    forwards[i],
                    strike,
                    t_exp,
                )

        return vols
