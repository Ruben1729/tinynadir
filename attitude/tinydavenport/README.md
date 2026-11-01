# tinydavenport

Part of [Attitude](../README.md) in [TinyNadir](../../README.md).

Davenport's q-method solves Wahba's weighted attitude-estimation problem using the maximum eigenpair of a symmetric 4×4 attitude-profile matrix. This is an alternative to TRIAD or QUEST, not another estimator correction to apply after them.

Include `tinydavenport/tinydavenport.h`; link `tinydavenport::tinydavenport` (with transitive `tinylinalg` and `tinyrotation` dependencies). Inputs are flattened arrays of corresponding reference/body direction vectors, with 2–64 observations. Vectors are normalized internally; nonnegative weights are rescaled without overflowing and normalized to sum to one. NULL weights select equal weights. Zero-weight pairs are ignored, though every supplied component must be finite.

The output quaternion uses scalar-first Hamilton convention and maps reference coordinates to body coordinates. The optional loss is `0.5*sum(normalized_weight*|body-C*reference|²)`. All positive-weight vectors must be nonzero, and the normalized maximum-eigenvalue gap must exceed 1e-12; otherwise the estimate is rejected as nonunique/numerically unobservable. This is not a mission-specific quality gate or covariance estimate. All vectors must use a common epoch and reference frame.

```c
tinydavenport_workspace work;
double q[4], loss;
tinydavenport_status status = tinydavenport_solve(count, reference_xyz,
    body_xyz, weights, q, &loss, &work);
```

Caller-owned workspace is **544 bytes** when `sizeof(double)==8`, independent of observation count. It must not overlap inputs/outputs. Outputs are preserved on error; workspace is scratch. The symmetric eigensolver uses at most 64 Jacobi sweeps; no heap, I/O, mutable globals, or unbounded loops occur. Bounded local arrays and libm stack are additional working stack, not included in the workspace figure.

## Verification

96 frozen independent **SciPy 1.18.1 Rotation.align_vectors** results use the Kabsch/SVD formulation rather than Davenport's eigenproblem. Seed 821780 covers 2, 3, 4, 8, 16, and 64 observations; varied weights; zero-weight observations; exact 180° rotations; and noiseless/noisy vectors. Tests compare quaternions modulo sign (3e-12 tolerance) and normalized loss (3e-13). They also cover equal/huge weights, collinear observations, all-zero/negative/NaN weights, invalid counts, NULL inputs, and error preservation.

Development generators/raw outputs/logs are ignored under `.validation/tinydavenport/`. Retained C fixtures/tests are offline and do not depend on Python, NumPy, SciPy, or AHRS.

References: [Davenport/Wahba mathematical description](https://ahrs.readthedocs.io/en/latest/filters/davenport.html), [SciPy align_vectors](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.align_vectors.html). Independently written implementation; root MIT license and module `NOTICE` apply.
