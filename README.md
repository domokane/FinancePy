[![unit test action](https://github.com/ru-corporate/FinancePy/actions/workflows/run-unit-tests.yml/badge.svg)](https://github.com/ru-corporate/FinancePy/actions/workflows/run-unit-tests.yml)

![alt text](./images/logo.jpg?raw=true)

# Latest News

2 November 2021: A new version 0.204 has just been released that can be installed via pip.

# DISCLAIMER

This software is distributed FREE & WITHOUT ANY WARRANTY. 

Report any bugs or suggestions here as an issue. 

# CONTRIBUTORS WANTED !

If you have a knowledge of Quantitative Finance and a reasonable knowledge of Python, then please consider contributing to this project. There are small tasks and big tasks to be done. Just look in the list of Issues and you may find something you can do. Before you begin, please comment in the issue thread in case someone else may be working on that issue. Or you can contact me directly at dominic.okane at edhec.edu. 

If you are a user and require some additional functionality, then please add it as an issue.

# Quick Start Guide

FinancePy can be installed from pip using the following command:

`pip install financepy`

To upgrade an existing installation type:

`pip install --upgrade financepy`

I have encountered problems using Anaconda3-2020.07 due to some Numba and LLVMLite problems. However Anaconda3-2020.02 works.

## Using FinancePy in a Jupyter Notebook

Once financepy has been installed, it is easy to get started.

Just download the project and examine the set of Jupyter Notebooks in the notebooks folder.

A pdf manual describing all of the functions can be found in the project directory.

## Overview

FinancePy is a python-based library that is currently in beta version. It covers the following functionality:

* Valuation and risk models for a wide range of equity, FX, interest rate and credit derivatives.

Although it is written entirely in Python, it can achieve speeds comparable to C++ by using Numba. As a result the user has both the ability to examine the underlying code and the ability to perform pricing and risk at speeds which compare to a library written in C++.

The target audience for this library includes:

* Students of finance and students of python
* Academics teaching finance or conducting research into finance
* Traders wishing to price or risk-manage a derivative.
* Quantitative analysts seeking to price or reverse engineer a price.
* Risk managers wishing to replicate and understand price sensitivity.
* Portfolio managers wishing to check prices or calculate risk measures.
* Fund managers wanting to value a portfolio or examine a trading strategy.

Users should have a good, but not advanced, understanding of Python. In terms of Python, the style of the library has been determined subject to the following criteria:

1. To make the code as simple as possible so that those with a basic Python fluency can understand and check the code.
2. To keep all the code in Python so users can look through the code to the lowest level.
3. To offset the performance impact of (2) by leveraging Numba to make the code as fast as possible without resorting to Cython.
4. To make the design product-based rather than model-based so someone wanting to price a specific product can easily find that without having to worry too much about the model – just use the default – unless they want to. For most products, a Monte-Carlo implementation has been provided both as a reference for testing and as a way to better understand how the product functions in terms of payments, their timings and conditions.
5. To make the library as complete as possible so a user can find all their required finance-related functionality in one place. This is better for the user as they only have to learn one interface.
6. To avoid complex designs. Limited inheritance unless it allows for significant code reuse. Some code duplication is OK, at least temporarily.
7. To have good documentation and easy-to-follow examples.
8. To make it easy for interested parties to contribute.

In many cases the valuations should be close to if not identical to those produced by financial systems such as Bloomberg. However for some products, larger value differences may arise due to differences in date generation and interpolation schemes. Over time it is hoped to reduce the size of such differences.

Important Note:

* IF YOU HAVE ANY PRICING OR RISK EXAMPLES YOU WOULD LIKE REPLICATED, SEND SCREENSHOTS OF ALL THE UNDERLYING DATA, MODEL DETAILS AND VALUATION.
* IF THERE IS A PRODUCT YOU WOULD LIKE TO HAVE ADDED, SEND ME THE REQUEST.
* IF THERE IS FUNCTIONALITY YOU WOULD LIKE ADDED, SEND ME A REQUEST.

## The Library Design

The underlying Python library is split into a number of major modules:

* Utils - These are utility functions used to assist you with modelling a security. These include dates (Date), calendars, schedule generation, some finance-related mathematics functions and some helper functions.
* Market - These are modules that capture the market information used to value a security. These include interest rate and credit curves, volatility surfaces and prices.
* Models - These are the low-level models used to value derivative securities ranging from Black-Scholes to complex stochastic volatility models.
* Products - These are the actual securities and range from Government bonds to Bermudan swaptions.

Any product valuation is the result of the following data design:

**VALUATION** = **PRODUCT** + **MODEL** + **MARKET**

The interface to each product has a value() function that will take a model and market to produce a price.

## Author

Dominic O'Kane. I am a Professor of Finance at the EDHEC Business School in Nice, France. I have 12 years of industry experience and 10 years of academic experience.

Contact me at dominic.okane at edhec.edu.

## Dependencies

FinancePy depends on Numpy, Numba, Scipy and basic python libraries such as os, sys and datetime.

## Changelog

See the changelog for a detailed history of changes.

## Contributions

Contributions are very welcome. There are a number of requirements:

* The code should be Pep8 compliant.
* Comments are required for every class and function and they should be a clear description.
* At least one broad test case and a set of unit tests must be provided for every function.
* Avoid very pythonic constructions. For example a loop is as good as a list comprehension. And with numba it can be faster. Readability is the priority.

## License

 GPL-3.0 License - See the license file in this folder for details.
