# Control

[TinyNadir](../README.md) / Control

Construct desired attitudes, calculate torque/dipole requests, allocate actuators, and evaluate rotational dynamics. Select only the pieces needed by the control mode; drivers and mode supervision remain application responsibilities.

## Assemble a control mode

| Mode | Possible chain | What is optional |
| --- | --- | --- |
| Detumbling | Calibrated magnetometer → `tinybdot` → application torquer driver | No orbit, IGRF, attitude estimate or wheels required |
| Wheel pointing | Attitude/rate input + `tinyguidance` → optional `tinyslew` → PD/PID → `tinywheels` | External attitude can replace onboard estimation |
| Magnetic pointing | Attitude/rate input → PD/PID → `tinymagnetic` | Requested torque parallel to the field is unavailable |
| Wheel unloading | Wheel momentum + field → `tinydesaturation` → coordinated wheel/torquer commands | Application owns mode coordination and counter-torque |
| Dynamics simulation | Applied/disturbance torques → `tinyrigidbody` derivatives → selected [integrator](../core/README.md) | Environment, sensors and actuator simulation are application choices |

[Guidance](tinyguidance/README.md) supplies desired motion; [attitude estimation](../attitude/README.md) supplies actual motion; a controller requests torque; allocation converts it to actuator commands. These are separate functions. Choose PD or PID as appropriate; actuator drivers and a supervisor remain application responsibilities.

Runnable examples: [minimal B-dot](../examples/detumble.c) and [wheel pointing](../examples/wheel_pointing.c). The wheel example accepts an attitude input without linking an orbit or environment model. Gains, actuator authority, sampling, sensor contamination, eclipse handling and termination criteria require mission-specific validation.

## Storage and model limits

Wheel allocation uses **312 bytes** of caller workspace; the rigid-body kernel uses **168 bytes**, plus bounded stack. Wheel/magnetic allocators project and clip; they do not solve an optimal constrained allocation or redistribute clipped wheel commands. Slews are rest-to-rest, acceleration-limited profiles without jerk/exclusion constraints. Rotational dynamics assumes fixed inertia and returns derivatives rather than integrating internally.

Independent fixtures use NumPy/SciPy matrix, rotation, filter and trajectory equivalents. Tests exercise signs, saturation, state preservation and relevant physical invariants. They validate the documented discrete laws and selected models rather than establishing closed-loop stability for arbitrary spacecraft.

## Algorithm glossary

| Algorithm | Purpose and supported profile | Inputs → outputs; connections |
| --- | --- | --- |
| [tinyguidance](tinyguidance/README.md) | Two-vector, target-line and nadir desired attitudes. | Same-frame geometry/orbit → desired quaternion and supported reference rate; separate from estimating actual attitude. |
| [tinyslew](tinyslew/README.md) | Shortest rest-to-rest attitude slew with triangular/trapezoidal speed. | Initial/target quaternions, rate/acceleration bounds, elapsed time → reference attitude/rate/acceleration; no jerk or exclusion constraints. |
| [tinybdot](tinybdot/README.md) | Filtered B-dot magnetic detumbling. | Consecutive calibrated body magnetic fields → bounded dipole; no orbit, IGRF or attitude solution required. |
| [tinypd](tinypd/README.md) | Quaternion PD with optional torque feedforward. | Desired/current attitude, current-body rate error, gains → requested body torque. |
| [tinypid](tinypid/README.md) | Quaternion PID with bounded integral, per-axis torque limits and conditional anti-windup. | Attitude/rate error and caller state → torque; alternative to PD. |
| [tinymagnetic](tinymagnetic/README.md) | Minimum-norm magnetic torque projection and component dipole clipping. | Requested torque and body field → dipole, achievable torque and residual; no instantaneous field-parallel authority. |
| [tinywheels](tinywheels/README.md) | Reaction-wheel torque allocation for 1–8 wheels. | Body torque, axes and torque/speed/enable limits → rotor torques/residual; pseudoinverse followed by clipping, without redistribution. |
| [tinydesaturation](tinydesaturation/README.md) | Magnetic wheel-momentum unloading toward an optional target. | Wheel momentum and body field → dipole/achievable external torque; application coordinates wheel counter-torque and modes. |
| [tinyrigidbody](tinyrigidbody/README.md) | Euler rigid-body angular acceleration and quaternion kinematics. | Fixed inertia, optional wheel momentum and applied body torque → derivatives; integrate separately. |
| [tinydisturbance](tinydisturbance/README.md) | Gravity-gradient, residual-dipole, lever-arm, drag-plate and SRP-plate torques/forces. | Supplied environment and body geometry → disturbance; optional simulation/analysis input. |

## Build an individual algorithm

From the collection root:

```sh
cmake -S control/tinybdot -B build/tinybdot -DCMAKE_BUILD_TYPE=Release
cmake --build build/tinybdot
ctest --test-dir build/tinybdot --output-on-failure
```

Keep the full checkout for shared CMake helpers and dependencies across groups. Link `tinybdot::tinybdot`; public include paths and C API names are unchanged. Building one algorithm includes its required dependencies, not every algorithm in this group. See the [collection build options](../README.md#build-and-test).
