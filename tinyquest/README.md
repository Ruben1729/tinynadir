# tinyquest

Heap-free C99 **QUEST (Quaternion Estimator)** for the weighted Wahba attitude problem. It constructs the 4×4 attitude profile, solves the maximum root of its characteristic quartic with bounded Newton iterations, and recovers the quaternion from the adjugate of the shifted profile. Full adjugate-column selection handles 180° rotations without the singularity of a scalar-first Rodrigues recovery. Three bounded Rayleigh/adjugate refinements reduce polynomial cancellation. This is not a wrapper around Davenport's eigensolver.

Include `tinyquest/tinyquest.h`; link `tinyquest::tinyquest`. Inputs are flattened arrays of 2–64 matching reference/body directions. Nonnegative weights are normalized after scaling; NULL means equal weights, and zero-weight pairs are ignored. All components must be finite; positive-weight directions must be nonzero. The scalar-first Hamilton quaternion maps reference coordinates to body coordinates, following `tinyrotation`. Loss, when requested, is normalized `0.5*sum(w*|body-C*reference|²)`.

```c
tinyquest_workspace work;
double q[4], loss;
tinyquest_status status = tinyquest_solve(count, reference_xyz, body_xyz,
                                        weights, q, &loss, &work);
```

Caller-owned workspace is **384 bytes** for an eight-byte double, regardless of observation count. Output quaternion/loss are preserved on error; workspace is scratch and must be separate from inputs/outputs. The maximum is 64 Newton steps plus three fixed refinements. Other bounded local arrays and math-library frames are additional stack. No heap, shared mutable globals, I/O, or external reference runtime is used.

An insufficient maximum-root derivative (<=1e-12), tiny adjugate recovery (column squared norm <1e-24), or unresolved optimal eigenvector returns a failure/degeneracy status. Nearly unobservable data need a mission-specific geometry/quality gate. Positive weights and finite vectors do not guarantee a usable attitude. Input vectors must correspond to the same epoch and common reference frame. A second weighted solver such as Davenport is an alternative module, not an independent measurement to fuse afterward.

## Verification

96 frozen **SciPy 1.18.1 Rotation.align_vectors** Kabsch/SVD solutions (NumPy 2.5.4, seed 821780) independently verify the rotation and normalized loss. They cover 2, 3, 4, 8, 16, and 64 observations; random directions/orientations; unequal/zero weights; noisy/noiseless cases; and exact 180° rotations about each axis. Tests also cover equal/huge weights, degeneracy, invalid weights/counts/NULLs, and failure-state preservation. Quaternion error is measured modulo sign with a 3e-12 limit; loss tolerance is 3e-13.

The development generator, raw independent results, and comparison log are ignored under `.validation/tinyquest/`. Only numeric fixtures and standalone C tests are committed; no Python dependency is required for normal verification.

References: [QUEST mathematical formulation](https://ahrs.readthedocs.io/en/latest/filters/quest.html), [Shuster and Oh, 1981](https://ntrs.nasa.gov/citations/19810033694), and [SciPy align_vectors](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.align_vectors.html). The implementation is independently written under the root MIT license; see `NOTICE`.
