# Market Volatility
These modules create a family of curve types related to the market volatility. There are three types of class:

1. Term structures of volatility i.e. volatility as a function of option expiry date.
2. Volatility curves which are smile/skews so store volatility as a function of option strike.
3. Volatility surfaces which hold volatility as a function of option expiry date AND option strike.

The classes are as follows:

### equity_vol_surface
Constructs an equity volatility surface that fits to a grid of market volatilities at a set of strikes and expiry dates. It implements the SVI parametric form for fitting and interpolating volatilities. It also provides plotting of the volatility curve and surfaces.

### FinFXVolSurface
FX volatility as a function of option expiry and strike. This class constructs the surface from the ATM volatility plus a choice of 10 and 25 delta strangles and risk reversals or both. This is done for multiple expiry dates. A number of curve fitting choices are possible including polynomial in delta and SABR.

## FinIborCapFloorVol
Libor cap/floor volatility as a function of option expiry (cap/floor start date). Takes in cap (flat) volatility and bootstraps the caplet volatility. This is assumed to be piecewise flat.

## FinIborCapFloorVolFn
Parametric function for storing the cap and caplet volatilities based on form proposed by Rebonato. 