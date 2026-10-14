# Contributing to FinancePy

Thank you for your interest in contributing to FinancePy!

FinancePy is an open-source Python library for the valuation and risk analysis of financial securities and derivatives. Its aim is to provide financial models that are accurate, readable, efficient, and accessible to practitioners, researchers, and students.

Contributions of all sizes are welcome, including bug fixes, new pricing models, improved numerical methods, additional tests, documentation, and examples.

## 1. Getting Started

Before contributing, please:

- Read the [README](https://github.com/domokane/FinancePy/blob/master/README.md) to understand the project's objectives and design philosophy.
- Explore the existing code, tests, and examples.
- Check the [open issues](https://github.com/domokane/FinancePy/issues) and [pull requests](https://github.com/domokane/FinancePy/pulls) to avoid duplicating work.
- For substantial changes or new models, consider opening an issue first to discuss the proposed approach.

Small bug fixes and documentation improvements can generally be submitted directly as pull requests.

## 2. Development Setup

Fork the repository on GitHub and clone your fork:

```bash
git clone https://github.com/YOUR_USERNAME/FinancePy.git
cd FinancePy
```

Create a virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate
```

On Windows, activate the environment using `.venv\Scripts\activate`.

Install the project and its dependencies:

```bash
python -m pip install -r requirements.txt
python -m pip install -e .
```

Create a branch for your changes:

```bash
git checkout -b feature/my-contribution
```

## 3. Coding Guidelines

FinancePy places a high priority on code that is easy to read, understand, and maintain.

Please follow these general principles:

- Follow PEP 8 conventions and the style of the surrounding code.
- Use descriptive names for variables, functions, and classes.
- Keep functions reasonably focused and avoid unnecessary complexity.
- Document the purpose, inputs, outputs, and assumptions of new public functionality.
- Reuse existing FinancePy utilities and conventions wherever possible.
- Avoid introducing additional dependencies unless there is a compelling reason.

**Readability is important.** FinancePy is intended not only for practical financial calculations but also for understanding how financial models work.

Where possible, mathematical implementations should be sufficiently transparent for a reader to relate the code to the underlying financial theory.

## 4. Financial and Numerical Correctness

Financial correctness is a fundamental requirement for contributions to FinancePy.

When introducing or modifying a pricing model, numerical algorithm, or financial convention, please:

- Explain the mathematical methodology and its assumptions.
- Provide references to relevant academic papers, textbooks, or established methodologies where appropriate.
- Validate results against analytical solutions, published benchmarks, or independently verified calculations.
- Test representative market conditions and important edge cases.
- Use numerical tolerances appropriate to the algorithm and problem.

Relevant edge cases may include zero volatility, very short maturities, extreme strikes, negative interest rates, or limiting parameter values, depending on the model.

Where appropriate, contributors should also test known financial relationships, such as put-call parity, convergence properties, or consistency between alternative valuation methods.

Differences arising from market conventions, interpolation methods, or modelling assumptions should be explained rather than hidden by overly generous test tolerances.

**Correctness and clarity should take priority over computational performance.**

## 5. Testing

New functionality should normally be accompanied by appropriate tests.

FinancePy contains a `tests/` directory with existing tests that can serve as examples.

There are two folders

- Regression tests are tests of overall functionality and changes in the output are monitored using differences.
- Unit tests are standard pytest compatible function tests.

Before submitting a pull request:

- Run the relevant existing tests.
- Add tests for new functionality and bug fixes.
- Ensure that existing behaviour has not been unintentionally changed.
- Include regression tests where a defect has been corrected.

Tests should be reproducible and should not depend on external market data services or other resources that may be unavailable to contributors.

For numerical calculations, avoid unnecessarily strict floating-point comparisons. Tolerances should reflect the accuracy expected from the underlying method.

Consult the repository's existing test scripts and development instructions for the appropriate commands.

## 6. Performance and Numba

FinancePy uses Numba to accelerate computationally intensive numerical routines.

Performance improvements are welcome, particularly for pricing algorithms, simulations, and numerical calculations that may be executed repeatedly.

However:

- Prefer a clear and correct implementation before optimizing it.
- Use Numba where it provides a measurable benefit and is consistent with the surrounding code.
- Avoid unnecessary allocations and repeated calculations in performance-critical routines.
- Ensure that optimizations do not introduce numerical instability or change results unexpectedly.
- Where possible, provide benchmarks when proposing significant performance improvements.

Remember that Numba compilation can introduce overhead on the first execution. Performance comparisons should distinguish compilation time from subsequent execution time.

## 7. Adding New Products or Models

Contributions introducing new financial products or pricing models are particularly welcome.

A new implementation should ideally include:

1. **A clear description** of the financial instrument or model.
2. **Mathematical references** explaining the pricing methodology.
3. **Documented assumptions**, conventions, and input parameters.
4. **Tests** against analytical solutions or reliable benchmarks.
5. **An example** demonstrating how the functionality can be used.

Please follow the existing organization of the library:

- `financepy/products/` — financial instruments and valuation interfaces.
- `financepy/models/` — mathematical and numerical models.
- `financepy/market/` — curves, volatility structures, and market data.
- `financepy/utils/` — dates, calendars, schedules, and mathematical utilities.

Where possible, new functionality should integrate with existing classes and interfaces rather than introducing parallel implementations of the same concepts.

Not every contribution needs to be extensive. A well-tested improvement to an existing model can be just as valuable as adding an entirely new product.

## 8. Documentation and Examples

Documentation improvements are always appreciated.

Useful contributions include:

- Improving docstrings and explanations.
- Adding worked valuation examples.
- Clarifying financial conventions and model assumptions.
- Correcting outdated documentation.
- Providing examples that illustrate how different FinancePy components work together.

Examples should be reasonably self-contained and easy for other users to reproduce.

## 9. Submitting a Pull Request

When your changes are ready:

1. Ensure that the relevant tests pass.
2. Review your changes for clarity and consistency.
3. Update documentation or examples where appropriate.
4. Push your branch to your fork.
5. Open a pull request against the FinancePy repository.

Please include a concise explanation of:

- What has changed and why.
- Any relevant mathematical or financial references.
- How the changes have been tested.
- Any limitations, assumptions, or compatibility considerations.

Keep pull requests focused where practical. Smaller, self-contained changes are generally easier to review and maintain.

Maintainers may request changes or additional tests before merging. This is part of ensuring the reliability and consistency of the library.

## 10. Reporting Bugs

If you discover a bug, please [open an issue](https://github.com/domokane/FinancePy/issues).

A useful bug report should include:

- The FinancePy and Python versions being used.
- A minimal example reproducing the problem.
- The actual result or error message.
- The expected result and why it is expected.
- Relevant financial conventions or model assumptions.

For numerical discrepancies, please include the input parameters and, where possible, an independent reference value.

## 11. Licensing

FinancePy is distributed under the GNU General Public License, version 3 or later (GPL-3.0-or-later).

By submitting a contribution, you agree that your contribution may be distributed under the project's license.

Please ensure that submitted code is your own work or that you have the necessary rights to contribute it under this license.

## 12. Final Remarks

FinancePy aims to make financial modelling more transparent and accessible while maintaining high standards of numerical accuracy and software quality.

Contributions from quantitative analysts, software developers, academics, students, and financial practitioners are all welcome.

Whether you are correcting a typo, improving an existing model, or implementing a new pricing method, your contribution can help make FinancePy more useful to the wider quantitative finance community.

**Thank you for helping improve FinancePy!**
