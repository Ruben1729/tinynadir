# tinycalibration

Apply a supplied three-axis affine sensor calibration, followed by alignment into spacecraft body axes:

`body = sensor_to_body * correction * (raw - bias)`

Include `tinycalibration/tinycalibration.h`; link `tinycalibration::tinycalibration`. Matrices are contiguous row-major arrays. The correction can represent scale, cross-axis, or magnetometer soft-iron corrections; the alignment must be a proper rotation. The module **applies** parameters: it does not estimate calibration from a data set, validate the physical appropriateness of the correction matrix, or own measurement timestamps. Units follow the supplied values/parameters.

Parameters and observations must be finite. Invalid alignment or arithmetic overflow returns a status. Output may alias the raw input and remains unchanged on failure. All storage is caller parameters/output plus bounded local arrays; no heap, I/O, mutable globals, or hardware dependencies are used.

## Verification

64 frozen affine transforms from **NumPy 2.5.4 matrix multiplication** with **SciPy 1.18.1 Rotation** alignment matrices (seed 331377) compare the complete chain. Maximum observed difference: 7.11e-15; tolerance is 5e-14 times max(1, expected magnitude). Reflection/nonrotation/NaN/NULL and output-preservation checks accompany the valid cases. These comparisons verify the specified affine operation, not sensor calibration quality.

Raw independent results/generators remain ignored under `.validation/tinycalibration/` and `.validation/adcs-reference/`; static C regression fixtures run offline. See `NOTICE` and `LICENSE`.
