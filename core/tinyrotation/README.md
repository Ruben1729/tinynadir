# tinyrotation

Part of [Core](../README.md) in [TinyNadir](../../README.md).

Heap-free C99 rotations using **scalar-first Hamilton quaternions** `[w,x,y,z]`. The numeric matrix `C(q)` maps a reference-frame vector into body-frame coordinates: `b=C(q)*r`. Quaternion multiplication satisfies `C(p*q)=C(p)*C(q)`. Matrices are row-major; angles use radians. This convention must be adapted when connecting a library that uses vector-first quaternions, JPL multiplication, or the inverse frame mapping.

Include `tinyrotation/tinyrotation.h`; link `tinyrotation::tinyrotation`. APIs cover quaternion normalization/composition/inversion, quaternion↔DCM, axis-angle and rotation-vector conversions, ZYX Euler angles, vector rotation, and shortest-arc SLERP. No heap, I/O, globals, or dynamically sized storage occurs; working arrays are bounded at nine doubles. Inputs/outputs may alias where the documented shapes permit it. Failed calls preserve output arrays.

```c
double q[4], body[3];
const double axis[3] = {0, 0, 1}, reference[3] = {1, 0, 0};
tinyrotation_from_axis_angle(axis, 1.5707963267948966, q);
tinyrotation_apply(q, reference, body); /* approximately {0,1,0} */
```

Input quaternions and axis directions are normalized using scaling to avoid overflow. A zero quaternion/axis is degenerate. DCM inputs must already be proper rotations within 1e-10; no hidden matrix repair is performed. Axis-angle/rotation-vector output selects the principal angle in [0,pi], with the x-axis chosen for a zero angle. Euler ZYX means `Rz(yaw)*Ry(pitch)*Rx(roll)`; at gimbal lock output roll is zero and yaw absorbs the unobservable combination. Other Euler sequences are outside this profile.

## Verification

**SciPy 1.18.1 Rotation and Slerp**, with NumPy 2.5.4, produced 96 frozen cases (seed 12914): small angles, principal axes, 90° and 180° rotations, nearly 180°, full turns, and random orientations. Static C tests compare DCMs, rotated vectors, rotation vectors, composition, and SLERP; they also check inverse/round trips, gimbal lock, antipodal quaternion interpolation, invalid matrices, NaNs, and huge quaternion normalization. Quaternion comparisons account for the equivalent signs `q` and `-q`.

Observed maximum DCM/SLERP differences and all pass/fail results are recorded in the development comparison log under ignored `.validation/tinyrotation/`. Tolerances are 3e-14 for DCM/vector/rotation-vector values and 1e-14 quaternion Euclidean error after sign selection; Euler round-trip tolerance is 3e-13. Reference generators, packages, and raw outputs remain ignored; normal tests require only C and libm.

Reference: [SciPy Rotation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.html). New code is independently written from rotation identities under the root MIT license; see `NOTICE`.
