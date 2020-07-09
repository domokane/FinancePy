#Libor Products
This folder contains a set of Libor-related products. More recently with the demise of Libor these are known as Ibor products. It includes:

## FinInterestRateFuture 
This is a class to handle interest rate futures contracts. This is an exchange-traded contract to receive or pay Libor on a specified future date. It can be used to build the Liboir term structure. 

##FinLiborCapFloor 
This is a contract to buy a sequence of calls or puts on Libor over a period at a strike agreed today.

##FinLiborDeposit 
This is the basic Libor instrument in which a party borrows an amount for a specified term and rate unsecured.

##FinLiborFRA
This is a class to manage Forward Rate Agreements (FRAs) in which one party agrees to lock in a forward Libor rate.

##FinLiborSwap
This is a contract to exchange fixed rate coupons for floating Libor rates. This class has functionality to value the swap contract and to calculate its risk.

##FinLiborSwaption
This is a contract to buy or sell an option on a swap. The model includes code that prices a payer or receiver swaption.

##FinOIS
This is a contract to exchange the daily compounded Overnight index swap rate for a fixed rate agreed at contract initiation.
