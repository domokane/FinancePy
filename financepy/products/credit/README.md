This folder contains a set of credit-related assets ranging from CDS to CDS options, to CDS indices, CDS index options and then to CDS tranches. They are as follows:
* FinCDS is a credit default swap contract. It includes schedule generation, contract valuation and risk-management functionality.
* FinCDSBasket is a credit default basket such as a first-to-default basket. The class includes valuation according to the Gaussian copula.
* FinCDSIndexOption is an option on an index of CDS such as CDX or iTraxx. A full valuation model is included.
* FinCDSOption is an option on a single CDS. The strike is expressed in spread terms and the option is European style. It is different from an option on a CDS index option. A suitable pricing model is provided which adjusts for the risk that the reference credit defaults before the option expiry date.
* FinCDSTranche is a synthetic CDO tranche. This is a financial derivative which takes a loss if the total loss on the portfolio exceeds a lower threshold K1 and which is wiped out if it exceeds a higher threshold K2. The value depends on the default correlation between the assets in the portfolio of credits. This also includes a valuation model based on the Gaussian copula model.

### FinCDSCurve
This is a curve that has been calibrated to fit the market term structure of CDS contracts given a recovery rate assumption and a FinIborSingleCurve discount curve. It also contains a IborCurve object for discounting. It has methods for fitting the curve and also for extracting survival probabilities.
