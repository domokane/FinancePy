# FinancePy

FinancePy is a python-based library that is currently in beta version. It covers the following functionality:

* Valuation and risk models for a wide range of equity, FX, interest rate and credit derivatives.

Although it is written entirely in Python, it can achieve speeds comparable to C++ by using Numba. As a result the user has both the ability to examine the underlying code and the ability to perform pricing and risk at speeds which compare to a library written in C++.

The target audience for this library includes:

* Students wishing to learn derivative pricing and Python.
* Professors wishing to teach derivative pricing and Python.
* Traders wishing to price or risk-manage a derivative.
* Quantitative analysts seeking to price or reverse engineer a price.
* Risk managers wishing to replicate and understand a price.
* Portfolio managers wishing to check prices or calculate risk measures
* Fund managers wanting to value a portfolio or examine a trading strategy
* Structurers or financial engineers seeking to examine the pricing of a derivative structure.

Users should have a good, but not advanced, understanding of Python. In terms of Python, the style of the library has been determined subject to the following criteria:

1. To make the code as simple as possible so that those with a basic Python fluency can understand and check the code.
2. To keep all the code in Python so users can look through the code to the lowest level.
3. To offset the performance impact of (2) by leveraging Numba to make the code as fast as possible without resorting to Cython.
4. To make the design product-based rather than model-based so someone wanting to price a specific product can easily find that without having to worry too much about the model – just use the default – unless they want to. For most products, a Monte-Carlo implementation has been provided both as a reference for testing and as a way to better understand how the product functions in terms of payments, their timings and conditions.
5. To make the library as complete as possible so a user can find all their required finance-related functionality in one place. This is better for the user as they only have to learn one interface.
6. To avoid complex designs. Some code duplication is OK, at least temporarily.
7. To have good documentation and easy-to-follow examples.
8. To make it easy for interested parties to contribute.

In many cases the valuations should be close to if not identical to those produced by financial systems such as Bloomberg. However for some products, larger value differences may arise due to differences in date generation and interpolation schemes. Over time it is hoped to reduce the size of such differences.

Important Note:

* IF YOU HAVE ANY PRICING OR RISK EXAMPLES YOU WOULD LIKE REPLICATED, SEND SCREENSHOTS OF ALL THE UNDERLYING DATA, MODEL DETAILS AND VALUATION.
* IF THERE IS A PRODUCT YOU WOULD LIKE TO HAVE ADDED, SEND ME THE REQUEST.
* IF THERE IS FUNCTIONALITY YOU WOULD LIKE ADDED, SEND ME A REQUEST.

## Quick Start Guide

FinancePy can be installed from pip using the command:

`pip install financepy`

To upgrade an existing installation type:

`pip install --upgrade financepy`

You can then access a set of Jupyter Notebooks on at Github

<https://github.com/domokane/FinancePy-Examples>

Note that the first import of the library may take over 10 seconds as Numba compiles all of the imported module functions with a Numba decorator - this makes certain python functions run at speeds comparable with the C programming language. These compiled functions are then cached locally. Thereafter FinancePy should load almost instantly.

Download a zipped version by clicking on the green button and selecting Download ZIP.

A pdf manual describing all of the functions can be found at the same repository.

## The Library Design

The underlying Python library is split into a number of major modules:

* Finutils - These are utility functions used to assist you with modelling a security. These include dates (FinDate), calendars, schedule generation, some finance-related mathematics functions and some helper functions.
* Market - These are modules that capture the market information used to value a security. These include interest rate and credit curves, volatility surfaces and prices.
* Models - These are the low-level models used to value derivative securities ranging from Black-Scholes to complex stochastic volatility models.
* Products - These are the actual securities and range from Government bonds to Bermudan swaptions.

Any product valuation is the result of the following data design:

**VALUATION** = **PRODUCT** + **MODEL** + **MARKET**

The interface to each product has a value() function that will take a model and market to produce a price.

## Author

Dominic O'Kane. I am a Professor of Finance at the EDHEC Business School in Nice, France. I have 12 years of industry experience and 10 years of academic experience.

## Dependencies

FinancePy depends on Numpy, Numba, Scipy and basic python libraries such as os, sys and datetime.

## Changelog

See the changelog for a detailed history of changes

## Contributions

Contributions are very welcome. There are a number of requirements:

* You should use CamelCase i.e. variables of the form optionPrice
* Comments are required for every class and function and they should be clear
* At least one test case must be provided for every function
* Follow the style of the code as currently written. This may change over time but please use the current style as your guide.

## Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## License

MIT
