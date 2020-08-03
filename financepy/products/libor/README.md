#Libor Products
This folder contains a set of Libor-related products. More recently with the demise of Libor these are known as Ibor products. It includes:

## FinInterestRateFuture 
This is a class to handle interest rate futures contracts. This is an exchange-traded contract 
to receive or pay Libor on a specified future date. It can be used to build the Liboir term structure. 

##FinLiborCapFloor 
This is a contract to buy a sequence of calls or puts on Libor over a period at a strike agreed today.

##FinLiborDeposit 
This is the basic Libor instrument in which a party borrows an amount for a specified term and rate unsecured.

##FinLiborFRA
This is a class to manage Forward Rate Agreements (FRAs) in which one party agrees to lock in a forward Libor rate.

##FinLiborSwap
This is a contract to exchange fixed rate coupons for floating Libor rates. This class has functionality to value the swap contract and to calculate its risk.

##FinLiborSwaption
This is a contract to buy or sell an option to enter into a swap to either pay or receive a fixed swap rate at a specific future expiry date. The model includes code that prices a payer or receiver swaption with the following models:
- Black's Model
- Shifted Black Model
- SABR
- Shifted SABR
- Hull-White Tree Model
- Black-Karasinski Tree Model
- Black-Derman-Toy Tree Model

##FinLiborBermudanSwaption
This is a contract to buy or sell an option to enter into a swap to either pay or receive a fixed swap rate at a specific future expiry date on specific coupon dates starting on a designated expiry date. The model includes code that prices a payer or receiver swaption with the following models:
- Hull-White Tree Model
- Black-Karasinski Tree Model
- Black-Derman-Toy Tree Model

It is also possible to price this using a Libor Market Model. However for the moment this must be done directly via the Monte-Carlo implementation of the LMM found in FinModelRatesLMM.

##FinOIS
This is a contract to exchange the daily compounded Overnight index swap rate for a fixed rate agreed at contract initiation.

## FinLiborCurve
This is a discount curve that is extracted by bootstrapping a set of Libor deposits, Libor FRAs and Libor swap prices. The internal representation of the curve are discount factors on each of the deposit, FRA and swap maturity dates. Between these dates, discount factors are interpolated according to a specified scheme - see below.
