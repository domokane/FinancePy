# Curves
These modules create a family of curve types related to the term structures of interest rates. There are two basic types of curve:

1. Best fit yield curves fitting to bond prices which are used for interpolation. A range of curve shapes from polynomials to B-Splines is available.
2. Discount curves that can be used to present value a future cash flow. These differ from best fits curves in that they exactly refit the prices of bonds or CDS. The different discount curves are created by calibrating to different instruments. They also differ in terms of the term structure shapes they can have. Different shapes have different impacts in terms of locality on risk management performed using these different curves. There is often a trade-off between smoothness and locality.

## Best Fit Bond Curves
The first category are BondYieldCurves.

### BondYieldCurve
This module describes a curve that is fitted to bond yields calculated from bond market prices supplied by the user. The curve is not guaranteed to fit all of the bond prices exactly and a least squares approach is used. A number of fitting forms are provided which consist of 

* Polynomial 
* Nelson-Siegel
* Nelson-Siegel-Svensson
* Cubic B-Splines

This fitted curve cannot be used for pricing as yields assume a flat term structure. It can be used for fitting and interpolating yields off a nicely constructed yield curve interpolation curve.

### FinCurveFitMethod
This module sets out a range of curve forms that can be fitted to the bond yields. These includes a number of parametric curves that can be used to fit yield curves. These include:
* Polynomials of any degree 
* Nelson-Siegel functional form. 
* Nelson-Siegel-Svensson functional form.
* B-Splines

## Discount Curves
These are curves which supply a discount factor that can be used to present-value future payments.

### FinDiscountCurve
This is a curve made from a Numpy array of times and discount factor values that represents a discount curve. It also requires a specific interpolation scheme. A function is also provided to return a survival probability so that this class can also be used to handle term structures of survival probabilities. Other curves inherit from this in order to share common functionality.

### FinDiscountCurveFlat
This is a class that takes in a single flat rate. 

### FinDiscountCurveNS
Implementation of the Nelson-Siegel curve parametrisation.

### FinDiscountCurveNSS
Implementation of the Nelson-Siegel-Svensson curve parametrisation.

### FinDiscountCurveZeros
This is a discount curve that is made from a vector of times and zero rates.


### FinInterpolate
This module contains the interpolation function used throughout the discount curves when a discount factor needs to be interpolated. There are three interpolation methods:

1. PIECEWISE LINEAR - This assumes that a discount factor at a time between two other known discount factors is obtained by linear interpolation. This approach does not guarantee any smoothness but is local. It does not guarantee positive forwards (assuming positive zero rates).
2. PIECEWISE LOG LINEAR - This assumes that the log of the discount factor is interpolated linearly. The log of a discount factor to time T is T x R(T) where R(T) is the zero rate. So this is not linear interpolation of R(T) but of T x R(T).
3. FLAT FORWARDS - This interpolation assumes that the forward rate is constant between discount factor points. It is not smooth but is highly local and also ensures positive forward rates if the zero rates are positive.
