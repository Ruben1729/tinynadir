# tinypolynomial

Part of [Core](../README.md) in [TinyNadir](../../README.md).

Bounded all-root solver for real polynomials of degree1–8, including complex roots. Coefficients descend in degree. Linear/quadratic cases use direct stable formulas; higher orders use a scaled simultaneous Aberth iteration with256 iteration limit and a backward-residual check. Caller owns16 complex workspace entries at degree8 (256bytes); bounded local storage is a nine-double coefficient array, eight output complexes, and scalars. No allocation, globals or I/O. Repeated/clustered roots can be ill-conditioned and may fail convergence or have large forward error despite a small residual. This is not a certified root-isolation solver; applications select physically meaningful Gauss/Lambert branches separately.

## Verification

96 seeded, separated real/conjugate-root polynomial cases of degree1–8 and multiple coefficient/root scales compare with independently computed NumPy companion-matrix roots. Root sets are matched without relying on conjugate order; scaled forward tolerance2e-8 in these conditioned cases. Zero roots, invalid coefficients/degrees and output atomicity are tested. Reference: https://numpy.org/doc/stable/reference/generated/numpy.roots.html .

Reference scripts, packages and raw outputs are in ignored `.validation/tinypolynomial/`. Regression fixtures are checked in and run offline without Python or reference libraries.
