# tinyrk4

Part of [Core](../README.md) in [TinyNadir](../../README.md).

Classical fixed-step fourth-order Runge–Kutta. A derivative callback makes the kernel usable for orbit, attitude, or other bounded ODE states (1–32 components). Six caller-provided state-length work arrays; no heap, globals or I/O. Four derivative calls per nonzero step. Backward steps are supported. Header specifies output atomicity and callback contract. Fixed steps have discretization error; this API does not estimate or enforce an error tolerance.

## Verification

Frozen independent SciPy DOP853 integration results compare exponential decay, harmonic motion, and nonlinear Lorenz dynamics over small positive and negative intervals. Per-case tolerances include the expected RK4 local truncation error, rather than claiming equality with a high-order solver. Raw source/version records are in the reference JSON; checked-in fixtures and tests print the sampled maximum scaled error. Reference: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html .

Reference scripts, packages and raw outputs are in ignored `.validation/tinyrk4/`. Regression fixtures are checked in and run offline without Python or reference libraries.
