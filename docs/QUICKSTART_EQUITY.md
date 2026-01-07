## Quickstart Guide to Equity Derivatives
We start with valuing a plain vanilla Equity Option and calculating its delta.

The easiest way to begin is to use a wildcard import from the utilities folder to pull in Dates and other required utility classes.
```
from financepy.utils import *
```
Then load in the equity derivative classes from the equity products folder. You will load classes you will not use but this is the simplest way to begin.
```
from financepy.products.equity import *
```
Now use the Date class to define the option expiry date for the 1st of June 2026
```
expiry_dt = Date(1, 6, 2026)
```
We then need to specify the option strike price
```
strike_price = 50.0
```
We finally create the option object. We are going to define an option with a payoff of a European call. A European call is a call option with a single expiry date.
```
call_option = EquityVanillaOption(expiry_dt, strike_price, OptionTypes.EUROPEAN_CALL)
```
There are alternative payoffs such as an AMERICAN_CALL. This can be exercised at any time before expiry.

We can print the option object to check that it is what we wanted
```
print(bond)
```
**Output:**
```
OBJECT TYPE: EquityVanillaOption
EXPIRY DATE: 01-JUN-2026
STRIKE PRICE: 50.0
OPTION TYPE VALUE: OptionTypes.EUROPEAN_CALL
NUMBER: 1.0
```
---
This data is exactly what you would find on the term sheet of this equity call option.

To perform a valuation we need additional market information starting with a stock price
```
stock_price = 55.0
```
We need a volatility
```
volatility = 0.20
```
We need a risk-free interest rate
```
interest_rate = 0.05
```
We need a dividend yield
```
dividend_yield = 0.02
```
and we need the valuation date. We choose the 1 Jan 2026 so that the option has exactly six months to expiry.
```
value_dt = Date(1, 1, 2026)
```
Now the valuation needs the interest rate in the form of a discount curve. We choose the simple flat curve at the interest rate provided
```
discount_curve = DiscountCurveFlat(value_dt, interest_rate, FrequencyTypes.CONTINUOUS)
```
We also need a curve for the dividend yield
```
dividend_curve = DiscountCurveFlat(value_dt, dividend_yield, FrequencyTypes.CONTINUOUS)
```
We need to specify a model. We use the standard Black-Scholes Model which takes in a single flat volatility
```
model = BlackScholes(volatility)
```
We can then call the model valuation of the call option as so
```
call_option.value(value_dt, stock_price, discount_curve, dividend_curve, model)
```
**Output:**
```
np.float64(2.8449436313648704)
```
The option is worth $2.845.

We can also obtain the option's delta
```
call_option.delta(value_dt, stock_price, discount_curve, dividend_curve, model)
```
**Output**
```
np.float64(0.559227735475412)
```
It is 0.559.

Finally, we can use vectorisation to plot the option price as a function of the stock price. We write
```
stock_prices = np.linspace(20,80,100)
```
This produces a vector of prices from $20 to $80 with 100 intervals. We then call the option pricing function again but passing in the vector of stock prices and store the output in variable $value$
```
value = call_option.value(value_dt, stock_prices, discount_curve, dividend_curve, model)
```
We can then plot this
```
plt.plot(stock_prices, value)
plt.xlabel("Stock Price")
plt.ylabel("Option Premium")
```
![alt text](image.png)

That's all. To see more, look at the notebooks in the notebooks folder.
