# Financial Validation Tests

This directory contains tests of FinancePy's financial correctness.

Validation tests compare FinancePy calculations against independent
reference results, such as:

- analytical solutions;
- published examples;
- market-standard reference calculations;
- independently implemented calculations;
- trusted external reference datasets.

These tests are separate from `tests/`, which primarily tests
implementation behavior, edge cases, and software correctness.

## Requirements

Every validation test should identify:

1. the financial quantity being validated;
2. the reference result;
3. the source or derivation of that result;
4. the numerical tolerance;
5. the reason for the chosen tolerance.

Expected values should not be produced by the FinancePy implementation
being tested.