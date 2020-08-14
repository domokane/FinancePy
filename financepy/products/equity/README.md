#Equity Products
This folder contains a set of Equity-related products. It includes:

## FinEquityVanillaOption
Handles simple European-style call and put options on a dividend paying stock with analytical and monte-carlo valuations.

## FinEquityAmericanOption
Handles America-style call and put options on a dividend paying stock with tree-based valuations.

## FinEquityAsianOption 
Handles call and put options where the payoff is determined by the average-stock price over some period before expiry.

## FinEquityBasketOption
Handles call and put options on a basket of assets, with an analytical and Monte-Carlo valuation according to Black-Scholes model.

## FinEquityCompoundOption
Handles options to choose to enter into a call or put option. Has an analytical valuation model for European style options and a tree model if either or both options are American style exercise.

## FinEquityDigitalOption
Handles European-style options to receive cash or nothing, or to receive the asset or nothing. Has an analytical valuation model for European style options.

## FinEquityFixedLookbackOption
Handles European-style options to receive the positive difference between the strike and the minimum (put) or maximum (call) of the stock price over the option life. 

## FinEquityFloatLookbackOption
Handles an equity option in which the strike of the option is not fixed but is set at expiry to equal the minimum stock price in the case of a call or the maximum stock price in the case of a put. In other words the buyer of the call gets to buy the asset at the lowest price over the period before expiry while the buyer of the put gets to sell the asset at the highest price before expiry. '''
    
## FinEquityBarrierOption
Handles an option which either knocks-in or knocks-out if a specified barrier is crossed from above or below, resulting in owning or not owning a call or a put option. There are eight variations which are all valued.

## FinEquityRainbowOption
TBD

## FinEquityVarianceSwap
TBD

Products that have not yet been implemented include:
* Power Options
* Ratchet Options
* Forward Start Options
* Log Options
