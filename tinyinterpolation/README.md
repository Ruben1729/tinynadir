# tinyinterpolation

Bounded scalar table interpolation: piecewise linear, up to eight-point Lagrange, and cubic Hermite with endpoint slopes. No extrapolation, mutable globals, heap or I/O. Header documents input limits; functions preserve outputs on failure. Use mission-relative seconds or otherwise well-scaled coordinates. Lagrange can be ill-conditioned on clustered points; users choose an appropriate interpolation window. Maximum local array is eight doubles. Link `tinyinterpolation::tinyinterpolation`.

## Verification

Frozen seeded cases compare NumPy interpolation and SciPy BarycentricInterpolator/CubicHermiteSpline, including derivatives. Relative-to-max(1, expected) limits: 2e-12 linear/Hermite, 2e-11 Lagrange. Source references: https://docs.scipy.org/doc/scipy/reference/interpolate.html and https://numpy.org/doc/stable/reference/generated/numpy.interp.html . Reference versions are recorded in tests/data/reference.h.

Reference scripts, packages and raw outputs are in ignored `.validation/tinyinterpolation/`. Regression fixtures are checked in and run offline without Python or reference libraries.
