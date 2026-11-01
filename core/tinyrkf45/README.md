# tinyrkf45

Part of [Core](../README.md) in [TinyNadir](../../README.md).

Adaptive classical Fehlberg 4(5) ODE integration, with fifth-order accepted states and componentwise infinity-norm error control. The caller provides 9*n doubles of workspace for 1–32 state components and explicit step/tolerance/attempt limits. No allocation, mutable globals or I/O. Positive/negative time spans supported; no event or dense-output interface. Outputs and statistics are preserved on failure; callback side effects remain caller-owned. Non-stiff systems only; an error tolerance is a local estimate, not a global accuracy guarantee.

## Verification

Independent SciPy DOP853 solutions of decay, oscillator and Lorenz systems are frozen into offline tests. Tested with absolute/relative local tolerances 1e-13, against reference-state tolerances 2e-9 scaled by max(1, expected). Tests exercise rejected/work-limited integrations, callback failures, invalid workspace and zero spans. Versions/cases and observed errors recorded in ignored reference outputs. Reference: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html .

Reference scripts, packages and raw outputs are in ignored `.validation/tinyrkf45/`. Regression fixtures are checked in and run offline without Python or reference libraries.
