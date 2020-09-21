This folder contains a suite of bond-related functionality across a set of files and classes. They are as follows:

* FinAnnuity is a stream of cashflows that is generated and can be priced.
* FinBond is a basic fixed coupon bond with all of the associated duration and convexity measures. It also includes some common spread measures such as the asset swap spread and the option adjusted spread.
* FinBondCallable is a bond that has an embedded call and put option. A number of rate models pricing functions have been included to allow such bonds to be priced and risk-managed.
* FinBondFuture is a bond future that has functionality around determination of the conversion factor and calculation of the invoice price and determination of the cheapest to deliver. 
* FinBondMarket is a database of country-specific bond market conventions that can be referenced. These include settlement days and accrued interest conventions.
* FinBondOption is a bond option class that includes a number of valuation models for pricing both European and American style bond options. Models for European options include a Lognormal Price, Hull-White (HW) and Black-Karasinski (BK). The HW valuation is fast as it uses Jamshidians decomposition trick. American options can also be priced using a HW and BK trinomial tree. The details are abstracted away making it easy to use.
* FinConvertibleBond enables the pricing and risk-management of convertible bonds. The model is a binomial tree implementation of Black-Scholes which allows for discrete dividends, embedded puts and calls, and a delayed start of the conversion option.
* FinFloatingNote enables the pricing and risk-management of a bond with floating rate coupons. Discount margin calculations are provided.
* FinMortgage generates the periodic cashflows for an interest-only and a repayment mortgage. 


## Conventions

* All interest rates are expressed as a fraction of 1. So 3% is 0.03.
* All notionals of bond positions are given in terms of a notional amount.
* All bond prices are based on a notional of 100.0.
* The face of a derivatives position is the size of the underlying position.

## Bond Curves
These modules create a family of curve types related to the term structures of interest rates. There are two basic types of curve:

1. Best fit yield curves fitting to bond prices which are used for interpolation. A range of curve shapes from polynomials to B-Splines is available.
2. Discount curves that can be used to present value a future cash flow. These differ from best fits curves in that they exactly refit the prices of bonds or CDS. The different discount curves are created by calibrating to different instruments. They also differ in terms of the term structure shapes they can have. Different shapes have different impacts in terms of locality on risk management performed using these different curves. There is often a trade-off between smoothness and locality.

### FinBondYieldCurve
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
