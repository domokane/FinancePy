# FinancePy Unit Test Suite

This folder contains the unti testing framework for FinancePy.

Unit tests are not intended to capture model output as this can change slightly over time.

Unit tests detect changes in static deterministic low level libraries both inside and outside FinancePy

Here are some examples of what should be tested here:
- Date generation
- Calendars and holiday adjustments
- Day count conventions
- Interpolation code
- Spline code

When it comes to models, unit tests should only examine edge cases where the output is sure and will never change:
- Zero volatility
- Put-Call parity
- Various limits of say in- and out-of-the-money
- Known limits like P=100 for a bond when y=c

Be creative but avoid cases where the values could drift. Normal model output should be tested in the regression tests


### Running the Unit Test Cases

You can run them manually by typing

    pytest --ignore legacy

This makes sure that it does not examine code in the legacy folder.

Over time the legacy folder items will be reincluded in the unit tests.

