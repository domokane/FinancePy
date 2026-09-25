# UNDER CONSTRUCTION 

I am currently updating examples to be one per product, model or utility 
They will initially be available as python scripts 
They are intended to be educational and also provide basic tests 
The notebooks will also be replaced by new product specific notebooks. 

# FinancePy Examples

This directory contains executable Python examples demonstrating how to use FinancePy.

The examples are intended to complement the FinancePy documentation and unit tests by showing how products, market data and valuation models are combined in practical applications.

FinancePy follows the general valuation framework:

```text
VALUATION = PRODUCT + MODEL + MARKET
```

The examples illustrate each of these components and how they interact.

## Getting Started

Install FinancePy using:

```bash
pip install financepy
```

or upgrade an existing installation with:

```bash
pip install --upgrade financepy
```

If you are working directly from the FinancePy source repository, run the examples from the repository environment so that they use the version of FinancePy you are developing.

For example:

```bash
python examples/equity/FinEquityVanillaOption.py
```

The exact filenames available depend on the current version of the repository.

## Directory Structure

Examples are organised by financial product or functionality.

Typical categories include:

```text
examples/
├── bonds/
├── credit/
├── equity/
├── fx/
├── market/
├── rates/
└── utils/
```

Each directory contains standalone Python programs illustrating functionality from the corresponding FinancePy modules.

## What the Examples Show

Examples may demonstrate:

* construction of financial products;
* construction of discount, dividend and other market curves;
* use of valuation models;
* calculation of prices and risk measures;
* comparison of alternative valuation methods;
* numerical convergence;
* sensitivity to model and market parameters;
* financial identities and limiting cases;
* plotting of prices and sensitivities.

Some examples deliberately perform additional calculations and consistency checks. These are useful for understanding the behaviour of a model, but they should not be confused with the formal FinancePy unit-test suite.

## Examples versus Tests

The `examples` directory is primarily educational and diagnostic.

An example should answer questions such as:

```text
How do I construct this product?

What market data does it require?

Which valuation models can I use?

How do I calculate its value?

How does its value change when an input changes?
```

Formal automated tests belong in the FinancePy test directories.

Tests should contain assertions and provide reproducible pass/fail checks. Examples may instead print tables, compare methods and generate plots to make model behaviour easier to inspect.

An example may therefore contain diagnostic calculations such as:

```python
print("MODEL VALUE :", model_value)
print("MC VALUE    :", mc_value)
print("DIFFERENCE  :", model_value - mc_value)
```

or parameter sweeps such as:

```python
for volatility in volatilities:
    value = option.value(...)
    print(volatility, value)
```

These calculations can subsequently motivate formal regression or unit tests.

## Typical Example

A simple FinancePy example generally follows four steps.

### 1. Define Dates and Contract Terms

```python
from financepy.utils.date import Date

value_dt = Date(1, 1, 2026)
expiry_dt = Date(1, 1, 2027)
```

### 2. Construct the Product

Import the required product class and specify its contractual terms.

```python
product = ...
```

### 3. Construct Market Data and Model

For example, a valuation might require:

```python
discount_curve = ...
dividend_curve = ...
model = ...
```

### 4. Value the Product

```python
value = product.value(
    value_dt,
    ...,
)

print("VALUE:", value)
```

More detailed examples may then investigate sensitivities, convergence or comparisons with alternative models.

## Parameter Analysis

For numerical models it is often useful to examine the behaviour of the valuation across a range of inputs.

For example:

```python
values = []

for volatility in volatilities:

    model = ...

    value = product.value(
        value_dt,
        ...,
        model,
    )

    values.append(value)
```

The results can then be printed or plotted.

Useful analyses include value against:

```text
stock price
strike
volatility
interest rate
dividend yield
credit spread
correlation
time to maturity
number of Monte Carlo paths
number of tree steps
```

depending on the product.

## Numerical Convergence

Examples involving Monte Carlo simulation, numerical integration or trees may include convergence studies.

For Monte Carlo:

```text
NUM PATHS          VALUE
------------------------
1,000
5,000
10,000
50,000
100,000
```

For lattice methods:

```text
NUM STEPS          VALUE
------------------------
25
50
100
250
500
1000
```

Such examples are useful both for understanding numerical accuracy and for detecting implementation problems.

## Financial Consistency Checks

Where possible, examples should also demonstrate known financial relationships.

These may include:

* put-call parity;
* option decomposition identities;
* limiting cases;
* monotonicity relationships;
* analytic versus numerical valuation;
* Monte Carlo versus closed-form valuation;
* convergence between discrete and continuous approximations.

These checks make examples useful both as demonstrations and as tools for investigating model behaviour.

## Running an Example

From the FinancePy repository root:

```bash
python examples/<category>/<example_file>.py
```

Alternatively, change to the relevant examples directory and run the file directly if that example supports direct execution.

Some examples generate plots and therefore require Matplotlib, which is included among FinancePy's dependencies.

## First-Run Performance

FinancePy uses Numba for a number of computationally intensive models.

The first execution of an example using a Numba-compiled function may therefore take longer while the function is compiled. Compiled functions are cached, so subsequent runs are normally much faster.

## Writing New Examples

New examples should aim to be:

* easy to read;
* self-contained;
* financially meaningful;
* reproducible;
* explicit about important assumptions;
* consistent with the FinancePy public API.

Prefer straightforward Python over unnecessarily compact constructions.

Where useful, an example should begin with a simple valuation before adding more detailed analysis.

A larger example might follow this structure:

```text
1. Product construction
2. Basic valuation
3. Alternative valuation methods
4. Sensitivity analysis
5. Numerical convergence
6. Financial consistency checks
7. Plots
```

Examples should use public FinancePy interfaces wherever possible rather than relying on internal implementation details.

## Example Naming

Use descriptive filenames that identify the FinancePy product or functionality being demonstrated.

For example:

```text
FinBond.py
FinBondConvertible.py
FinEquityAsianOption.py
FinEquityRainbowOption.py
FinEquitySwap.py
```

Where a single file contains several demonstrations of the same product, keep them together when this makes it easier for users to explore that product.

## Contributions

Contributions of new examples are welcome.

A useful contribution should demonstrate a realistic application of FinancePy and should be understandable to someone familiar with Python and the relevant area of finance.

If an example exposes unexpected behaviour or a possible bug, consider adding a corresponding unit test as well.

## Further Documentation

The main FinancePy repository contains additional documentation and quick-start material covering bonds, equity derivatives, rate derivatives, credit derivatives and FX derivatives.

For the full library, installation instructions and source code, see the main FinancePy repository.

## Disclaimer

FinancePy is provided for educational, research and analytical use. Users should independently verify results before relying on them for financial decisions.
