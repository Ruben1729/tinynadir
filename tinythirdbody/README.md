# tinythirdbody

Differential point-mass third-body gravity in a common inertial frame, using SI units. The indirect origin acceleration is included. Small spacecraft/body-distance ratios use a cancellation-resistant logarithmic/expm1 form. Works with any externally supplied body ephemeris and nonnegative gravitational parameter. No heap, global mutable state, or I/O; local vector storage is 48 bytes plus scalars.

## Verification

Independent 90-digit mpmath point-mass gravity differences are frozen for 100 seeded Sun/Moon-like and stress geometries. Scaled acceleration tolerance is 2e-13 relative to max(1e-12, component). Invalid/collision/zero-baseline cases also tested. Reference primitive: https://mpmath.org/doc/current/ .

Reference scripts, packages and raw outputs are in ignored `.validation/tinythirdbody/`. Regression fixtures are checked in and run offline without Python or reference libraries.
