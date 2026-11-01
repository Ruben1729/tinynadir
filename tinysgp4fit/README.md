# tinysgp4fit

Heap-free C99 local differential correction of conventional SGP4 mean elements from Cartesian **TEME** observations. Include `tinysgp4fit/tinysgp4fit.h`, link `tinysgp4fit::tinysgp4fit`; dependencies are `tinybatchls` and `tinytle`/`tinysgp4`.

This implements an actual weighted nonlinear fit, not an element formatter or a seed-copy adapter. It evaluates SGP4 for candidate elements, minimizes weighted residuals, and reports fitted elements, local covariance, normalized RMS, iteration count, and model-evaluation count. Metadata/epoch are fixed. Results can initialize a propagator using `tinytle_initialize_elements` without producing or parsing TLE text.

```c
#include "tinysgp4fit/tinysgp4fit.h"
/* seed is tinytle_elements; samples are tinysgp4fit_observation[] in TEME. */
tinysgp4fit_options options;
tinysgp4fit_default_options(&options);
options.fit_bstar = 1;
/* Allocate caller-owned workspace/result and populate samples before calling:
   status = tinysgp4fit_fit(&seed, samples, count, &options, &workspace, &result); */
```

The fitted parameters are inclination, RAAN, eccentricity, argument of perigee, mean anomaly, mean motion, and optionally B*. Epoch, catalog identifier, ephemeris type, ndot/nddot, station biases, maneuvers and XP elements are not fitted. Inclination/RAAN/perigee/anomaly and covariance angular units are **degrees**; mean motion is **revolutions/day**, B* **inverse Earth radii**. These are format-native GP inputs, not osculating Keplerian elements. Output covariance is packed `parameter_count × parameter_count` in that order; a fixed B* has no covariance row.

Samples contain position metres, coordinate velocity metres/second, component standard deviations, and minutes relative to the fixed seed epoch. Times must be within ±10080 minutes. Select positions-only or positions plus velocities through `use_velocity`. Up to 10 epochs are supported (at least 2 with velocities, 3 positions-only). Sigmas must be positive for included components; they represent independent measurement errors. Transform/whiten differently framed or correlated measurements before fitting.

The solver requires a nearby valid seed and informative observations. Exactly circular/equatorial classical-element fits are not identifiable in this parameterization. Seeds/trials with eccentricity outside `(1e-12,1-1e-12)` or inclination within 1e-6 degree of 0/180 are rejected, and unresolved rank is reported as `UNOBSERVABLE`. B* fitting also rejects a zero-span arc. Weak geometry, insufficient span, or model mismatch can cause nonconvergence or unreliable local covariance; successful local convergence is not a global orbit-determination guarantee.

Internally eccentricity uses a logistic parameter and mean motion a logarithm, keeping them physical. Angular model evaluations wrap to 360 degrees. Parameter preconditioning, bounded central differences, and a model-specific relative cost tolerance 1e-8 account for finite-precision SGP4 objective evaluations. Limits are at most 100 outer iterations and `1+(2*n+1)*max_iterations+2*n` model evaluations, each propagating the selected epochs. Covariance is the inverse weighted local normal matrix with known supplied variances, transformed back from internal parameters; it is not residual-scaled.

No allocation, file I/O, or mutable globals are used. The caller supplies scratch storage (**10000 bytes on the tested 64-bit ABI with eight-byte double**) containing the least-squares workspace, scalar observations/sigmas, satellite and output state. SGP4/solver working stack and libm runtime are additional. Scratch is modified on failure; the result object is preserved. The workspace satellite is scratch, not a promised final fitted propagator—initialize explicitly from `result.elements` for later use.

Validation: **seven independent SciPy 1.18.1 LM fits using python-sgp4 2.27 C++ predictions**, and **28 frozen holdout states**. Coverage includes six/seven fitted parameters, positions-only/position+velocity fits, a low-eccentricity near-polar LEO and a resonant deep-space case. Noisy synthetic datasets have 8 epochs each; seed residual cost exceeds 1e9 for the position+velocity fits and drops to approximately 3–4, demonstrating actual correction.

On AppleClang 21 with strict C99 and ASan/UBSan, maximum holdout component differences were **3.62e-5 m** and **1.64e-8 m/s**; normalized covariance-component difference **5.39e-6**; weighted-cost difference **9.72e-9**. Static regression limits are 1 mm,1e-6 m/s,1e-4 normalized covariance, and 1e-5 weighted cost. Tests also verify metadata/fixed-B* preservation, residual reduction, changed seed parameters, invalid weights, zero-span B* arcs, circular-element degeneracy, iteration exhaustion, and result preservation. These are implementation-comparison differences, not orbit-tracking accuracy guarantees.

Raw reference outputs/generator stay ignored under `.validation/tinysgp4fit/`; the Python environment and all oracle dependencies are ignored too. Normal builds and regression tests use only the checked-in C fixtures. See LICENSE/NOTICE.
