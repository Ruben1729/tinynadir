# tinyguidance

Part of [Control](../README.md) in [TinyNadir](../../README.md).

Geometric spacecraft pointing references in C99, without heap allocation. Include `tinyguidance/tinyguidance.h`; link `tinyguidance::tinyguidance`.

`tinyguidance_point` builds a desired reference-to-body quaternion from two requested reference directions and two body pointing axes using TRIAD. Primary alignment has priority and the secondary direction fixes roll. This supports Sun, inertial, and other vector-defined pointing without owning a Sun/IGRF/orbit model. `tinyguidance_target` constructs its primary direction by subtracting spacecraft position from target position, with all positions/vectors in the same frame.

`tinyguidance_nadir` receives inertial position [m] and velocity [m/s]. Desired body z points toward Earth, x follows the velocity projected perpendicular to nadir, and y completes a right-handed frame (opposite orbit angular momentum). Its physical planar-reference body rate is `[0,-|r cross v|/|r|²,0]`. It does not infer orbit-plane precession from unavailable acceleration/jerk inputs. Do not supply Earth-fixed velocity as though it were inertial.

These routines construct desired attitude; they neither estimate actual attitude nor create a constrained slew. Use `tinyslew` and a controller separately. Zero/parallel geometry, invalid values, and arithmetic failure preserve outputs. TRIAD's 1e-10 separation-sine floor applies; stricter application geometry thresholds may be needed. Scratch is bounded local vectors/matrices, with no persistent state, I/O, or global mutation.

## Verification

64 independent **NumPy 2.5.4 orbital-frame constructions** and **SciPy 1.18.1 matrix-to-quaternion** results (seed 231217) compare nadir orientations and physical rates. Maximum quaternion difference is 3.72e-16 modulo sign (3e-13 tolerance); rates match within 3e-15 rad/s. Static tests also confirm Earth direction maps onto desired body +z, planar-rate sign, target direction, radial/zero geometry, NaNs, and output preservation. Basic TRIAD comparisons are covered by that module's separate AHRS fixture suite.

Raw generator/JSON/logs remain ignored under `.validation/tinyguidance/` and `.validation/adcs-reference/`; offline C tests retain numeric expectations only. See `NOTICE` and `LICENSE`.
