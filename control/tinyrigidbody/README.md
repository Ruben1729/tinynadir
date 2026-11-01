# tinyrigidbody

Part of [Control](../README.md) in [TinyNadir](../../README.md).

C99 Euler rigid-body rotational dynamics and quaternion kinematics. Include `tinyrigidbody/tinyrigidbody.h`; link `tinyrigidbody::tinyrigidbody`.

Given body physical angular rate omega, symmetric positive-definite body inertia J, optional body-coordinate wheel momentum h, and applied torque tau, compute `J*omega_dot=tau-omega cross (J*omega+h)`. The torque must already include wheel motor reaction; subtracting wheel momentum-rate again would count that reaction twice. The scalar-first reference-to-body quaternion derivative follows `tinygyro`: `qdot=-0.5*[0,omega]*q`.

All units are SI. The inertia tensor is checked through Cholesky before solving. The module returns **derivatives**, not an integrated trajectory; select an integrator separately. A fixed body inertia and the supplied instantaneous wheel momentum are assumed. Time-varying inertia and flexible dynamics are outside this profile.

Caller workspace contains 21 doubles (168 bytes with eight-byte double); bounded pivot/vector/quaternion temporaries and libm calls add stack. Inputs/output buffers/workspace must obey the documented separation; qdot and omega_dot outputs are distinct and preserved on error. There is no heap, persistent state, I/O, or mutable global storage.

## Verification

64 frozen **NumPy 2.5.4 linear-solve/cross-product** Euler dynamics cases include nonspherical/off-diagonal inertias, applied torques, and wheel momentum. Quaternion derivatives are independently evaluated by central differences of **SciPy 1.18.1 Rotation**. Angular-acceleration tolerance is 3e-13 rad/s² (maximum sampled difference 6.94e-17); quaternion derivative tolerance is 3e-10.

A separate asymmetric torque-free ten-second trajectory is frozen from **SciPy DOP853 solve_ivp**, rtol 2.3e-14/atol 1e-15. A test-only fixed-step RK4 run checks those states, conservation of kinetic energy, and inertial angular momentum, including the physical quaternion sign. Invalid inertia and error preservation are tested. These references verify the equations/integration interface rather than qualifying a spacecraft dynamics model.

Raw oracles, trajectory JSON, and logs remain ignored under `.validation/tinyrigidbody/` and `.validation/adcs-reference/`; retained tests are static C. See `NOTICE` and `LICENSE`.
