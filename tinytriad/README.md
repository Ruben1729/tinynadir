# tinytriad

TRIAD instantaneous attitude determination from two reference vectors and their corresponding body observations. The first vector is matched exactly; the second determines the remaining rotation. With noisy measurements this is a deterministic construction, not a weighted optimum. Use QUEST/Davenport when a weighted Wahba solution is needed.

Include `tinytriad/tinytriad.h`; link `tinytriad::tinytriad`, which supplies `tinylinalg` and `tinyrotation` transitively. Vectors must refer to the same two physical directions at a consistent epoch. Magnitudes are ignored. The returned row-major DCM maps reference to body; quaternions are scalar-first Hamilton in the `tinyrotation` convention.

```c
double q[4];
tinytriad_status status = tinytriad_solve(sun_reference, magnetic_reference,
                                         sun_body, magnetic_body, q);
```

The input vectors need not be unit length; normalization is scaled for large/small magnitudes. A zero vector or either pair with absolute sine of separation at most 1e-10 is rejected. This threshold bounds a geometric singularity and is not a sensor-quality guarantee. Nearly parallel Sun/magnetic directions should be gated by the application at a stricter mission-specific threshold when necessary. Missing Sun data cannot be replaced with zeroes.

All storage is bounded stack-local scratch and caller outputs; no heap, I/O, mutable globals, or persistent state exists. Outputs are preserved on failure. Vector references/models and sensor drivers are supplied by the application; TRIAD itself has no orbit, Sun, or IGRF dependency.

## Verification

128 frozen matrices from **AHRS 0.4.0 TRIAD** (NumPy 2.5.4, SciPy 1.18.1, seed 88741) cover uniformly randomized orientations and vector geometries, both noiseless pairs and independent vector perturbations. C tests compare matrices and equivalent quaternion matrices to 3e-14 per component and check first-vector alignment, degenerate/near-parallel geometry, invalid/NULL inputs, output preservation, and extreme input scales.

The ignored `.validation/tinytriad/` directory holds the development generator, raw reference JSON, and comparison log; only frozen numerical fixtures and offline C tests are retained in Git. No AHRS implementation is copied into the runtime.

References: [AHRS TRIAD mathematical description](https://ahrs.readthedocs.io/en/latest/filters/triad.html), [Shuster and Oh, 1981](https://ntrs.nasa.gov/citations/19810033694). The implementation is independently written; see the root MIT license and module `NOTICE`.
