# FinancePy Source-Code Examples

This directory collects the existing FinancePy source-code examples and groups them by financial area. The calculations were preserved from the original example files; this reorganisation does not invent new API usage.

## Running an example

From the repository root, run an example with Python, for example:

```console
python examples/code/bonds/example_bond.py
```

## Categories

- [`bonds/`](./bonds/) — Bonds and inflation-linked products (20 scripts)
- [`rates/`](./rates/) — Interest rates, Ibor/OIS, swaps and swaptions (18 scripts)
- [`equity/`](./equity/) — Equity products and volatility (18 scripts)
- [`fx/`](./fx/) — Foreign-exchange products and volatility (12 scripts)
- [`credit/`](./credit/) — CDS and credit products (10 scripts)
- [`market/`](./market/) — Discount curves and market-curve examples (9 scripts)
- [`models/`](./models/) — Pricing models and process simulations (24 scripts)
- [`utils/`](./utils/) — Dates, calendars, schedules, maths and other utilities (12 scripts)

## Notes

- These files originated as FinancePy executable examples/regression-style scripts, so some are larger than tutorial snippets.
- Run examples with FinancePy installed and with the project dependencies installed.
