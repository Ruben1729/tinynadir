# tinydesaturation

Part of [Control](../README.md) in [TinyNadir](../../README.md).

Magnetic reaction-wheel momentum unloading, separate from B-dot body detumbling. Include `tinydesaturation/tinydesaturation.h`; link `tinydesaturation::tinydesaturation` (which supplies `tinymagnetic`). Inputs are wheel angular momentum and optional desired momentum in body coordinates [N m s], a nonnegative gain [1/s], body magnetic field [T], and magnetorquer dipole limits [A m²].

The external torque request is `-gain*(momentum-target)`; a NULL target means zero. Magnetic allocation projects onto available field authority and clips dipole components. The function returns dipole, actual achievable external unloading torque, and an optional unachievable torque residual. The application must coordinate reaction-wheel counter-torque with its attitude controller and mode/power limits. Magnetic torque cannot remove instantaneous momentum parallel to the field; time evolution of the field matters.

Storage is bounded local vectors and caller outputs; no persistent state, heap, I/O, global state, or drivers exist. Output arrays must be distinct. Nonfinite values, negative gain, invalid magnetic allocation, or overflow return status without modifying outputs.

## Verification

64 frozen **NumPy 2.5.4 SVD magnetic pseudoinverse** reference calculations (seed 949221) compare limited dipole, actual torque, and residual. Maximum sampled dipole difference is 1.73e-16 A m²; limits are 3e-13 A m² and 3e-17 N m. Tests also check unloading direction, absent target, degenerate zero field, and preservation of outputs. This verifies the allocation/law, not complete closed-loop unloading behavior.

Development oracles/JSON/logs remain ignored under `.validation/tinydesaturation/` and `.validation/adcs-reference/`; frozen C regression tests are offline. See `NOTICE` and `LICENSE`.
