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
