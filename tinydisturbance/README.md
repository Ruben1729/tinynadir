# tinydisturbance

Heap-free C99 environmental torque models, with vectors in one consistent body frame and all quantities in SI. Include `tinydisturbance/tinydisturbance.h`; link `tinydisturbance::tinydisturbance`.

Available profiles:

- Gravity gradient: `3*mu/|r|³ * cross(unit(r),J*unit(r))`, with positive-definite symmetric inertia.
- Residual magnetic dipole: `cross(dipole,B)`.
- Generic force lever arm: `cross(center_of_pressure_arm,force)`.
- One-sided projected-area drag plate: force opposite atmosphere-relative spacecraft velocity, effective area `area*max(0,dot(normal,unit(velocity)))`.
- One-sided optical radiation-pressure plate with supplied specular/diffuse reflectivities. The supplied Sun direction points spacecraft-to-Sun; force follows incident photons in the opposite direction. Back faces receive zero force. The reflectivities are nonnegative and sum to at most one.

Atmospheric density, Sun pressure/direction, geometry/coefficients, and spacecraft state are inputs. This module does not fetch a density/ephemeris model or apply self-shadowing, accommodation models, complex geometry, or electromagnetic disturbance models. Multiple plates can be evaluated and summed by a caller. It is distinct from a translational propagator or controller.

All arrays are fixed local/caller storage; Cholesky validation uses bounded 3×3 scratch. Output force and torque buffers must be distinct; invalid inputs, degenerate directions, and overflow preserve outputs. No heap, I/O, mutable globals, or persistent state occurs.

## Verification

64 frozen **NumPy 2.5.4** vector/matrix evaluations (seed 738112) compare every model, including arbitrary body frames, nonspherical inertia, exposed/back faces, and lever arms. Tolerance is 3e-12 times max(1e-15, expected component magnitude); maximum sampled solar-torque error is 1.32e-23 N m. Separate tests check drag dissipation, zero field/torque geometries, perfectly reflecting/absorbing plate limits, invalid optical parameters, and error preservation. These are independently evaluated mathematical equivalents, not a claim of measured environmental-model accuracy.

Raw oracle scripts/results/logs stay ignored under `.validation/tinydisturbance/` and `.validation/adcs-reference/`; frozen C tests run offline. See `NOTICE` and `LICENSE`.
