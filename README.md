# FinancePy

FinancePy is a library of native Python functions which covers the following functionality:

* Valuation and risk models for a wide range of equity, FX, interest rate and credit derivatives.
* Portfolio asset allocation using Markovitz and other methods.

It is written entirely in Python. As a result the user has the ability to examine the underlying code and its logic. The target audience for this library is intended to include:

* Students wishing to learn derivative pricing and Python.
* Professors wishing to teach derivative pricing and Python.
* Traders wishing to price or risk-manage a derivative.
* Quantitative analysts seeking to price or reverse engineer a price.
* Risk managers wishing to replicate and understand a price.
* Portfolio managers wishing to check prices or calculate risk measures
* Fund managers wanting to value a portfolio or examine a trading strategy
* Structurers or financial engineers seeking to examine the pricing of a derivative structure.

Users are expected to have a good, but not advanced, understanding of Python. Up until now my main focus has been on financial derivatives. In general my approach has been:

1. To make the code as simple as possible so that those with a basic Python fluency can understand and check the code.
2. To keep all the code in Python so users can look through the code to the lowest level.
3. To offset the performance impact of (2) by leveraging Numba to make the code as fast as possible without resorting to Cython.
4. To make the design product-based rather than model-based so someone wanting to price a specific exotic option can easily find that without having to worry too much about the model – just use the default – unless they want to.
5. To make the library as complete as possible so a user can find all their required finance-related functionality in one place. This is better for the user as they only have to learn one interface.
6. To avoid complex designs. Some code duplication is OK, at least temporarily.
7. To have good documentation and easy-to-follow examples.
8. To make it easy for interested parties to contribute.

In many cases the valuations should be close to if not identical to those produced by financial systems such as Bloomberg. However for some products, larger value differences may arise due to differences in date generation and interpolation schemes. Over time it is hoped to reduce the size of such differences.

IF YOU HAVE ANY EXAMPLES YOU WOULD LIKE REPLICATED, SEND SCREENSHOTS OF ALL THE UNDERLYING DATA AND MODEL DETAILS

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

FinancePy can be installed using pip (see instructions below). A set of template Jupyter notebooks can be found under the github repository called FinancePy-Examples. The link is as follows:

https://github.com/domokane/FinancePy-Examples

A pdf manual describing all of the functions can be found at the same repository.

## Help Needed

The current version of the code is a beta. If you have any questions or issues then please send them to me. Contact me via the github page.

## Author

Dominic O'Kane. I am a Professor of Finance at the EDHEC Business School in Nice, France. I have 12 years of industry experience and 10 years of academic experience.

## Installation

FinancePy can be installed from pip using the command:

pip install financepy

To upgrade an existing installation type:

pip install --upgrade financepy

## Dependencies

FinancePy depends on Numpy, Numba and Scipy.

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
