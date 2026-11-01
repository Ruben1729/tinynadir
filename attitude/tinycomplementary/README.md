# tinycomplementary

Part of [Attitude](../README.md) in [TinyNadir](../../README.md).

An explicit nonlinear complementary attitude observer with gyro bias estimation, independent of any terrestrial gravity assumption. Include `tinycomplementary/tinycomplementary.h`; link `tinycomplementary::tinycomplementary`.

For normalized measured body directions and their modeled common-reference directions, compute the weighted error `sum(weight*(measured cross predicted))`. Bias evolves as `bias_dot=-ki*error`, is component-limited, and the corrected physical body rate is `gyro-new_bias+kp*error`. Scalar-first Hamilton reference-to-body quaternion kinematics use `qdot=-0.5*[0,corrected_rate]*q`. The discrete update is first-order Euler plus normalization, matching the explicit-observer profile; it is not the exact held-rate exponential from `tinygyro_step`.

The caller owns seven-double quaternion/bias state (56 bytes for binary64). Positive dt, finite nonnegative gains/bias limit, and finite data are required. Up to 64 reference/body pairs are accepted; weights are nonnegative, rescaled, and normalized. Zero count or all-zero weights gives gyro-only prediction. Active direction pairs must be nonzero. Inputs must have consistent frame/time; the caller gates missing, disturbed, or stale measurements. A single direction cannot establish full instantaneous attitude. An orbital accelerometer is not a substitute for a Sun/nadir reference.

Initialize through `tinycomplementary_initialize`. Failure preserves state. All scratch is bounded stack-local vectors; no heap, I/O, mutable globals, covariance, or reference-model dependency is present. Gains and integration sample rate require mission validation; first-order propagation is intended for suitably small rate-times-step increments.

## Verification

128 frozen samples from **AHRS 0.4.0 Mahony.updateIMU** (seed 551104, NumPy 2.5.4) independently compare quaternion and bias evolution. The oracle's body-to-reference quaternion is conjugated at the boundary. Its unit-z reference is treated mathematically as a vector observation; the CubeSat implementation does not assume that accelerometers measure gravity. Maximum sampled quaternion difference is recorded by the comparison log; tolerance is 3e-13 modulo sign, with the same absolute bias tolerance.

Additional tests demonstrate unobservable single-vector rotation, correction toward a two-vector known attitude with zero gyro, bias limiting, gyro-only prediction, invalid counts/timing, and atomic failure. The zero-gyro correction case is deliberately supported even though AHRS's wrapper returns early for exactly zero gyro. Raw oracle packages/scripts/JSON/logs remain ignored under `.validation/tinycomplementary/` and `.validation/adcs-reference/`; static C tests run offline.

Reference: [Mahony explicit complementary observer](https://ahrs.readthedocs.io/en/latest/filters/mahony.html). The equations are implemented independently; see `NOTICE` and `LICENSE`.
