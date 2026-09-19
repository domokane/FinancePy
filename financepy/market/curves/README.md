# Discount Curves
These modules create a family of discount curve types related to the term structures of interest rates. Discount curves that can be used to present value a future cash flow. These differ from best fits curves in that they exactly refit the prices of bonds or CDS. The different discount curves are created by calibrating to different instruments. They also differ in terms of the term structure shapes they can have. Different shapes have different impacts in terms of locality on risk management performed using these different curves. There is often a trade-off between smoothness and locality. These are curves which supply a discount factor that can be used to present-value future payments.

### DiscountCurve
This is a curve made from a Numpy array of times and discount factor values that represents a discount curve. It also requires a specific interpolation scheme. A function is also provided to return a survival probability so that this class can also be used to handle term structures of survival probabilities. Other curves inherit from this in order to share common functionality.

### DiscountCurveFlat
This is a class that takes in a single flat rate.

### DiscountCurveNS
Implementation of the Nelson-Siegel curve parametrisation.

### DiscountCurveNSS
Implementation of the Nelson-Siegel-Svensson curve parametrisation.

### DiscountCurveZeros
This is a discount curve that is made from a vector of times and zero rates.


### Interpolate
This module contains the interpolation function used throughout the discount curves when a discount factor needs to be interpolated.
