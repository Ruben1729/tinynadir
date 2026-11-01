# tinypd

Part of [Control](../README.md) in [TinyNadir](../../README.md).

Quaternion proportional/derivative attitude feedback producing **requested body torque**, in N m. Include `tinypd/tinypd.h`; link `tinypd::tinypd`.

Both attitudes use scalar-first Hamilton reference-to-body mapping. Define `q_error=q_current*inverse(q_desired)` and select its shortest representative by the sign of its scalar. The law is `tau=2*kp*sign(q_error.w)*q_error.xyz-kd*rate_error+feedforward`. Rate error is measured minus desired physical angular rate, **already expressed in current body axes**. Desired-rate frame conversion belongs to the caller. Gains are finite nonnegative per-axis values, not tuned automatically; optional feedforward is a caller-specified body torque.

For a small positive physical z displacement with a stationary target, the quaternion z error and commanded restoring torque are negative. The quaternion signs q and -q produce the same command. Exact 180-degree shortest-path direction is intrinsically ambiguous; application guidance should avoid relying on an arbitrary sign there.

No actuator saturation/allocation is hidden in this function. Use wheel/magnetic allocation separately. No heap, persistent state, I/O, or mutable globals exist. Invalid quaternions/inputs and overflow return status and preserve outputs.

## Verification

64 frozen **SciPy 1.18.1 Rotation** error-quaternion compositions and **NumPy 2.5.4** feedback evaluations (seed 331377) compare torque to 3e-14 N m. Maximum sampled difference is 2.22e-16 N m. Physical restoring direction, quaternion sign equivalence, negative gains/NULLs, and failure preservation are checked. This verifies the control law implementation; closed-loop stability requires mission inertia, timing, gains, and actuator validation.

Raw reference scripts/results/logs remain ignored under `.validation/tinypd/` and `.validation/adcs-reference/`. Retained C tests are offline. See `NOTICE` and `LICENSE`.
