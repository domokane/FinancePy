This folder contains a suite of bond-related functionality across a set of files and classes. They are as follows:

* Bond is a basic fixed coupon bond with all of the associated duration and convexity measures. It also includes some common spread measures such as the asset swap spread and the option adjusted spread.
* BondZero is a zero coupon bond. This is a bond issued at a deep discount that matures at par. Accrued interest is calculated by interpolating the price growth.
* BondAnnuity is a stream of cash flows that is generated and can be priced.
* BondCallable is a bond that has embedded call and put options. A number of rate models pricing functions have been included to allow such bonds to be priced and risk-managed.
* BondConvertible enables the pricing and risk-management of convertible bonds. The model is a binomial tree implementation of Black-Scholes which allows for discrete dividends, embedded puts and calls, and a delayed start of the conversion option.
* BondFRN enables the pricing and risk-management of a bond with floating rate coupons. Discount margin calculations are provided.
* BondFuture is a bond future that has functionality around determination of the conversion factor and calculation of the invoice price and determination of the cheapest to deliver. 
* BondMarket is a database of country-specific bond market conventions that can be referenced. These include settlement days and accrued interest conventions.
* BondMortgage generates the periodic cash flows for an interest-only and a repayment mortgage. 
* BondOption is a bond option class that includes a number of valuation models for pricing both European and American style bond options. Models for European options include a Lognormal Price, Hull-White (HW) and Black-Karasinski (BK). The HW valuation is fast as it uses Jamshidians decomposition trick. American options can also be priced using a HW and BK trinomial tree. The details are abstracted away making it easy to use.
* BondPortfolio is a portfolio of bonds.
* Yield Curve is a class to handle bond yield curves. It uses a variety of shapes to best-fit a set of bond yields.
* Zero curve is a class to perform an exact fit to a set of provided bonds using a piece-wise flat zero rate.

## Conventions
* All interest rates are expressed as a fraction of 1. So 3% is 0.03.
* All notionals of bond positions are given in terms of a notional amount.
* All bond prices are based on a notional of 100.0.
* The face of a derivatives position is the size of the underlying position.

## Bond Curves
These modules create a family of curve types related to the term structures of interest rates. These are best fit yield curves fitting to bond prices which are used for interpolation. A range of curve shapes from polynomials to B-Splines is available. This module describes a curve that is fitted to bond yields calculated from bond market prices supplied by the user. The curve is not guaranteed to fit all of the bond prices exactly and a least squares approach is used. A number of fitting forms are provided which consist of 

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
