# Market Volatility
These modules create a family of curve types related to the market volatility. There are three types of class:

1. Term structures of volatility i.e. volatility as a function of option expiry date.
2. Volatility curves which are smile/skews so store volatility as a function of option strike.
3. Volatility surfaces which hold volatility as a function of option expiry date AND option strike.

The classes are as follows:

### FinEquityVolCurve
Equity volatility as a function of option strike. This is usually a skew shape.

### FinFXVolSurface
FX volatility as a function of option expiry and strike. This class constructs the surface from the ATM volatility and 25 delta strangles and risk reversals and does so for multiple expiry dates.

## FinIborCapFloorVol
Libor cap/floor volatility as a function of option expiry (cap/floor start date). Takes in cap (flat) volatility and boostraps the caplet volatility. This is assumed to be piecewise flat.

## FinIborCapFloorVolFn
Parametric function for storing the cap and caplet volatilities based on form proposed by Rebonato. 