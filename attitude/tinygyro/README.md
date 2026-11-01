# tinygyro

Part of [Attitude](../README.md) in [TinyNadir](../../README.md).

Held-body-rate quaternion propagation and continuous quaternion derivatives, in C99 without heap allocation. Include `tinygyro/tinygyro.h` and link `tinygyro::tinygyro`.

The scalar-first Hamilton quaternion maps reference coordinates into body coordinates. Gyro inputs are **physical body angular velocity relative to that reference**, expressed in body axes, in rad/s. Thus `qdot=-0.5*[0,omega]*q`: a positive physical z-axis rotation from identity produces a negative quaternion z component. Bias is subtracted in the same axes. This sign differs from a quaternion that maps body to reference.

`tinygyro_step` applies the exact exponential for a corrected body rate held constant over a finite nonnegative step in seconds. It is not a coning integrator or a substitute for sampling varying rates. `tinygyro_derivative` provides derivatives for a separate integrator. Inputs are finite; zero quaternions, invalid pointers/times, and arithmetic overflow return status. Outputs may alias quaternion input and are preserved on failure. Scratch is bounded stack-local storage; no persistent state, I/O, globals, or allocation occurs.

## Verification

64 frozen **SciPy 1.18.1 Rotation** compositions (NumPy 2.5.4, seed 331377) compare exact held-rate propagation, with a maximum quaternion difference of 5.15e-16 modulo sign. Central differences of SciPy trajectories independently check the continuous derivative to 3e-10. Tests include the physical positive-z sign, bias correction, invalid input/time, and error preservation. Propagation tolerance is 3e-14. Normal C tests require no reference packages.

Development reference scripts, JSON, and logs live in ignored `.validation/tinygyro/` and `.validation/adcs-reference/`; only frozen numerical fixtures are committed. [SciPy Rotation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.html) provides the independent rotation equivalent. See `NOTICE` and `LICENSE`.
