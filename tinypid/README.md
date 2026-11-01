# tinypid

Quaternion PID attitude feedback with bounded integral state, conditional anti-windup, and per-axis torque clipping. Include `tinypid/tinypid.h`; link `tinypid::tinypid`. The module transitively uses `tinypd`/`tinyrotation` conventions: scalar-first reference-to-body Hamilton quaternions; physical measured-minus-desired rate error in current body axes.

The proportional/derivative/feedforward term is `tinypd`'s law. Integral state accumulates `2*sign(q_error.w)*q_error.xyz*dt` by rectangular integration, clips to supplied integral limits, and freezes an axis's proposed increment when it would push an already saturated requested torque further into saturation. Final torque is componentwise limited. This profile does not incorporate actuator-allocation residual feedback or coupled-axis optimal anti-windup.

Zero-initialize the caller's three-double state. All gains and limits are finite/nonnegative; dt is finite/nonnegative. Changed integral limits clamp the old state on a successful call. A saturation bitmask is optional. Invalid inputs/overflow preserve state and outputs. No heap, mutable globals, I/O, or hardware drivers occur. Controller gains, torque bounds, desired-rate transformation, and actuator integration remain application responsibilities.

## Verification

128 frozen sequential cases use **SciPy 1.18.1 quaternion compositions** and independently evaluated **NumPy 2.5.4 discrete anti-windup equations** (seed 949221). Both integral history and limited torque are compared to 3e-14; maximum sampled torque difference is 1.04e-17 N m. Tests include extended saturation with no integral windup, bounds, invalid dt, and atomic state/output failure.

The reference calculations verify this specified discrete law, not equivalence to every controller called PID or a spacecraft stability guarantee. Raw oracle data/generators/logs stay ignored under `.validation/tinypid/` and `.validation/adcs-reference/`; static C fixtures run offline. See `NOTICE` and `LICENSE`.
