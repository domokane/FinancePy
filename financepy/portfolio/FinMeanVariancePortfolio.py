import matplotlib.pyplot as plt


class MeanVariancePortfolio():

    def __# Let's define some functions we will use later

    def portfolioVariance(weights):
        weights = np.array(weights)
        var = np.dot(weights.T, np.dot(assetCovariance * 252, weights))
        return var
    
    def portfolioVolatility(weights):
        return np.sqrt(np.dot(weights.T, np.dot(assetCovariance * 252, weights)))
    
    def portfolioReturn(weights):
        return np.sum(assetReturns * weights) * 252
    
    def portfolioSharpeRatio(weights):
        return (portfolioReturn(weights) - rfr) / portfolioVolatility(weights)


# In[23]:


# Risk free rate 
rfr = 0.015


# In[24]:


# We generate different portfolio weightings and store their mean return and volatility
portfolioReturns = []
portfolioVolatilies = []


# ## Mean-Variance Portfolio Optimisation

# In[25]:


# We wish to use SciPy's optimization function
import scipy.optimize as sco


# ### Calculating the Efficient Frontier

# We need to loop over returns and find the portfolio that minimises the portfolio volatility. 

# In[26]:


# This means that we need to include a constraint that enforces the return given by variable ret
# With the constraint to have the weights sum to 1 we have two constraints which we write as follows
cons = ({'type': 'eq', 'fun': lambda x:  portfolioReturn(x) - tret},
        {'type': 'eq', 'fun': lambda x:  np.sum(x) - 1})


# In[27]:


minRet = min(assetReturns*252)
maxRet = max(assetReturns*252)
trets = np.linspace(minRet, maxRet, 50)
tvols = []
initialWeights = np.ones(numAssets)
bnds = tuple((0, 1) for x in initialWeights)
    
for tret in trets:
    cons = ({'type': 'eq', 'fun': lambda x:  portfolioReturn(x) - tret},
            {'type': 'eq', 'fun': lambda x:  np.sum(x) - 1})

    res = sco.minimize(portfolioVolatility, initialWeights, method='SLSQP', bounds=bnds, constraints=cons)
    frontierWeights = res['x']
    frontierRet = portfolioReturn(frontierWeights)
    frontierVol = portfolioVolatility(frontierWeights)    
    tvols.append(res['fun'])


# In[21]:


import matplotlib.pyplot as plt


# In[30]:


plt.figure(figsize=(16, 8))
plt.scatter(tvols, trets, c=(trets-rfr) / tvols, marker='o')
plt.grid(True)
plt.xlabel('Expected volatility')
plt.ylabel('Expected return')
plt.colorbar(label='Sharpe ratio')


# In[ ]:




