# FinancePy Source-Code Examples

This directory collects the existing FinancePy source-code examples and groups them by financial area. The calculations were preserved from the original example files; this reorganisation does not invent new API usage.

## Running an example

From the repository root, run an example with Python, for example:

```console
python examples/code/bonds/test_bond.py
```

On Windows, you can also double-click an example `.py` file when `.py` files are associated with Python. When launched directly from Windows Explorer, the example keeps its console open and waits for **Enter** before closing. Normal command-line, IDE, CI, and non-Windows execution is not paused.

## File naming

Example scripts use Python snake_case and the `test_<class_or_feature>.py` convention, for example `test_bond.py`, `test_equity_vanilla_option.py`, and `test_black_scholes.py`.

## Categories

- [`bonds/`](./bonds/) — Bonds and inflation-linked products (20 scripts)
- [`rates/`](./rates/) — Interest rates, Ibor/OIS, swaps and swaptions (18 scripts)
- [`equity/`](./equity/) — Equity products and volatility (18 scripts)
- [`fx/`](./fx/) — Foreign-exchange products and volatility (12 scripts)
- [`credit/`](./credit/) — CDS and credit products (10 scripts)
- [`market/`](./market/) — Discount curves and market-curve examples (9 scripts)
- [`models/`](./models/) — Pricing models and process simulations (24 scripts)
- [`utils/`](./utils/) — Dates, calendars, schedules, maths and other utilities (12 scripts)

## Support files

`fin_test_cases.py`, `helpers.py`, and `add_fp_to_path.py` remain in this directory because the existing examples depend on them. `double_click_pause.py` implements the Windows Explorer console behaviour. Supporting regression output and data files are kept alongside the category that uses them.

## Notes

- These files originated as FinancePy executable examples/regression-style scripts, so some are larger than tutorial snippets.
- Run examples from a FinancePy development environment with the project dependencies installed.
- The category folders are organisational only; the underlying financial calculations have not been rewritten.
