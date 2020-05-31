# FinancePy

FinancePy is a library of native Python functions which covers the following functionality:

* Valuation and risk of a wide range of equity, FX, interest rate and credit derivatives.
* Valuation models for a range of bonds including callable and puttable bonds.
* Portfolio risk measures for portfolios of the securities above.
* Optimal Portfolio asset allocation using Markovitz and other methods.
* Time series analysis of financial data using econometric techniques.

The target audience for this library is intended to include:

* Students, professor, or other academics seeking to teach derivative or asset pricing pricing via a library which permits the ability to dig into the underlying functions.
* Traders wishing to price or risk manage a derivative.
* Quantitative analysts
* Risk managers both on the buy and sell side
* Portfolio managers wishing to check prices or calculate risk measures
* Fund managers wanting to value a portfolio or examine a trading strategy
* Structurers or financial engineers seeking to examine the pricing of a derivative structure.

Users would be expected to have a good, but not advanced, understanding of Python, financial derivatives and some mathematics. 

Up until now my main focus has been on financial derivatives. Subsequent versions will also include asset selection, portfolio-level risk management, regulatory calculations and market analysis tools. In general my objectives have been:

1. To make the code as simple as possible so that students and those with a basic Python fluency can understand and check the code.
2. To keep all the code in Python so users can look through the code to the lowest level.
3. To offset the performance impact of (2) by leveraging Numba to make the code as fast as possible without resorting to Cython.
4. To make the design product-based rather than model-based so someone wanting to price a specific exotic option can easily find that without having to worry too much about the model – just use the default – unless they want to.
5. To make the library as complete as possible so a user can find all their required finance-related functionality in one place. This is better for the user as they only have to learn one interface.
6. To avoid complex designs as I do not want to make it too hard for unskilled Python programmers to use the library.
7. To have good documentation and easy-to-follow examples.
8. To make it easy for interested parties to contribute.

In many cases the valuations should be close to if not identical to those produced by financial systems such as Bloomberg. However for some products, larger value differences may arise due to differences in date generation and interpolation schemes. Over time I expect to reduce the size of such differences.

## The Library Design
The underlying Python library is split into a number of major modules:

* Finutils - These are utility functions used to assist you with modelling a security. These include dates (FinDate), calendars, schedule generation, some finance-related mathematics functions and some helper functions.
* Market - These are modules that capture the market information used to value a security. These include interest rate and credit curves, volatility surfaces and prices.
* Models - These are the low-level models used to value derivative securities ranging from Black-Scholes to complex stochastic volatility models. 
* Products - These are the actual securities and range from Government bonds to Bermudan swaptions.

Any price is the result of a PRODUCT + MODEL + MARKET. The interface to each product has a value() function that will take a model and market to produce a price.

There are also two other folders which are currently fairly empty: They are:
* Portfolio - This will be where portfolio allocation will go,
* Risk - This is for portfolio risk analysis


## How to Use the Library

FinancePy can be installed using pip (see instructions below). I have provided a range of template Jupyter notebooks under the github repository called FinancePy-Examples. The link is as follows:

https://github.com/domokane/FinancePy-Examples

A pdf description of functions can be found at the same repository.

## Help Needed

The current version of the code is very much a beta. Hence there is no guarantee on its exactness. If you have any questions or issues then please send them to me as a matter of urgency and I will do my best to investigate as quickly as possible.

## Author

My name is Dr. Dominic O'Kane. I am a finance professor at the EDHEC Business School in Nice, France.

## Installation

FinancePy can be installed from pip

pip install financepy

or to upgrade

pip install --upgrade financepy

## Dependencies

FinancePy depends on Numpy and Numba and Scipy.

## Changelog

See the changelog for a detailed history of changes

## Contributions

Contributions are very welcome. There are a number of requirements:

* You should use CamelCase i.e. variables of the form optionPrice
* Comments are required for every class and function and they should be clear
* At least one test case must be provided for every function
* Follow the style of the code as currently written. This may change over time but please use the current style as your guide.

## License

MIT
