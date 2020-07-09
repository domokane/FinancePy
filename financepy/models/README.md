# Models
## Overview
This folder contains a range of models used in the various derivative pricing models implemented in the product folder. These include credit models for valuing portfolio credit products such as CDS Tranches, Monte-Carlo based models of stochastics processes used to value equity, FX and interest rate derivatives, and some generic implementations of models such as a tree-based Hull White model. Because the models are useful across a range of products, it is better to factor them out of the product/asset class categorisation as it avoids any unnecessary duplication. 

In addition we seek to make the interface to these models rely only on fast types such as floats and integers and Numpy arrays.

These modules hold all of the models used by FinancePy across asset classes. 

The general philosophy is to separate where possible product and models so that these models have as little product knowledge as possible. 

Also, Numba is used extensively, resulting in code speedups of between x 10 and x 100.

# Generic Arbitrage-Free Models
There are the following arbitrage-free models:
* FinModelBlack is Black's model for pricing forward starting contracts (in the forward measure) assuming the forward is lognormally distributed.
* FinModelBlackShifted is Black's model for pricing forward starting contracts (in the forward measure) assuming the forward plus a shift is lognormally distributed. CHECK
* FinModelBachelier prices options assuming the underlying evolves according to a Gaussian (normal) process.
* FinSABR Model is a stochastic volatility model for forward values with a closed form approximate solution for the implied volatility. It is widely used for pricing European style interest rate options, specifically caps and floors and also swaptions.
* FinSABRShifted Model is a stochastic volatility model for forward value with a closed form approximate solution for the implied volatility. It is widely used for pricing European style interest rate options, specifically caps and floors and also swaptions.

The following asset-specific models have been implemented:

# Equity Models
* FinHestonModel 
* FinHestonModelProcess
* FinProcessSimulator

# Interest Rate Models

### Equilibrium Rate Models
There are two main short rate models.
* FinCIRRateModel is a short rate model where the randomness component is proportional to the square root of the short rate. This model implementation is not arbitrage-free across the term structure.
* FinVasicekRateModel is a short rate model that assumes mean-reversion and normal volatility. It has a closed form solution for bond prices. It does not have the flexibility to fit a term structure of interest rates. For that you need to use the more flexible Hull-White model.

### Arbitrage Free Rate Models
* FinBlackKaraskinskiRateModel is a short rate model in which the log of the short rate follows a mean-reverting normal process. It refits the interest rate term structure. It is implemented as a trinomial tree and allows valuation of European and American-style rate-based options.
* FinHullWhiteRateModel is a short rate model in which the short rate follows a mean-reverting normal process. It fits the interest rate term structure. It is implemented as a trinomial tree and allows valuation of European and American-style rate-based options. It also implements Jamshidian's decomposition of the bond option for European options.

# Credit Models
* FinGaussianCopula1FModel is a Gaussian copula one-factor model. This class includes functions that calculate the portfolio loss distribution. This is numerical but deterministic.
* FinGaussianCopulaLHPModel is a Gaussian copula one-factor model in the limit that the number of credits tends to infinity. This is an asymptotic analytical solution.
* FinGaussianCopulaModel is a Gaussian copula model which is multifactor model. It has a Monte-Carlo implementation.
* FinLossDbnBuilder calculates the loss distribution.
* FinMertonCreditModel is a model of the firm as proposed by Merton (1974).

# FX Models
