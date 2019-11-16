# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 17:54:03 2019

@author: Dominic O'Kane
"""


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def test_FinMeanVarPortfolio():

    from financepy.portfolio.FinMeanVarPortfolio import FinMeanVarPortfolio

    m = [0.05, 0.1, 0.12, 0.18]
    C = [[0.0064, 0.00408, 0.00192, 0],
         [0.00408, 0.0289, 0.0204, 0.0119],
         [0.00192, 0.0204, 0.0576, 0.0336],
         [0, 0.0119, 0.0336, 0.1225]]

    portfolio = FinMeanVarPortfolio(["A", "B", "C", "D"])
    portfolio.setReturnMeanCovar(m, C)

    retIn = 0.09
    wts = portfolio.optimize(retIn, 0.0, 1.0)
    retOut = portfolio.portfolioReturn(wts)
    volOut = portfolio.portfolioVolatility(wts)
    print("Return In:", retIn, "Out:", retOut, "Vol:", volOut)

    rets, vols = portfolio.efficientFrontier(0.06, 0.18, 20)

    rfr = 0.015
    if 1 == 1:
        plt.plot(vols, rets, '-o')
        plt.grid(True)
        plt.xlabel('Expected volatility')
        plt.ylabel('Expected return')

    dowJones30Df = pd.read_pickle('.\\data\\DowJonesPricesFromYahoo.pkl')
    dowJones30Df.dropna(inplace=True)

    returns = dowJones30Df.pct_change(periods=1)
    returns.dropna(inplace=True)
    returns.head()

    # We need to make a list of the tickers that we want
    assetList = ['AAPL','MSFT','IBM','GS','V','WMT','VZ','CSCO']
    
    numAssets = len(assetList)
    
    # Now we use this to select columns from the dataframe
    newReturns = returns[assetList]
    assetReturns = newReturns.mean()
    assetCovariance = newReturns.cov()
    assetCorrelations = newReturns.corr()
    
    rfr = 0.015

    portfolio = FinMeanVarPortfolio(assetList)
    portfolio.setReturns(returns)


#test_FinMeanVarPortfolio()
