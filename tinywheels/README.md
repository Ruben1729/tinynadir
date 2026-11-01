# tinywheels

Bounded reaction-wheel torque allocation for 1–8 wheels. Include `tinywheels/tinywheels.h`; link `tinywheels::tinywheels`. Each wheel axis is a unit body-frame triple in a contiguous per-wheel array. Positive motor torque accelerates the rotor along its axis; spacecraft reaction torque is its negative.

This profile computes the minimum-norm unconstrained allocation using an eigen-based pseudoinverse of the 3×3 axis Gram matrix, clips each motor torque, and blocks commands that would further exceed a positive/negative speed limit. Disabled wheels and wheels with zero torque/speed authority are excluded. Unreachable/limited torque is returned as `requested - achieved`. There is **no iterative redistribution after clipping** and no claim of an optimal constrained solution. Rank-deficient layouts can still return their reachable projection.

Workspace is caller-owned (312 bytes for eight-byte doubles), independent of wheel count. Bounded local storage includes eight torque commands and vectors. Outputs/workspace must use distinct buffers as documented; invalid axes/limits/counts, no available actuator, or arithmetic failure preserves outputs. Axis normalization tolerance is 1e-10; Gram eigenvalues at most 1e-12 of the largest are discarded. No heap, I/O, globals, or wheel drivers are included.

## Verification

96 frozen **NumPy 2.5.4 SVD pseudoinverse** cases (seed 949221) cover every wheel count 1–8, arbitrary axes, disabled wheels, rank-deficient arrays, torque clipping, and speed saturation. Both allocated motor torque and resulting residual are compared to 3e-12 N m; maximum sampled motor-torque difference is 1.37e-17 N m. Tests explicitly verify rotor/body reaction sign, overspeed authority, unavailable actuators, invalid counts, and output preservation.

Raw generators/reference JSON/logs remain ignored under `.validation/tinywheels/` and `.validation/adcs-reference/`. Retained C tests are independent of NumPy and use static numerical expectations. See `NOTICE` and `LICENSE`.
