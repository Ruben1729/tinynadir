# Mission

[TinyNadir](../README.md) / Mission

Navigation, orbit fitting, uncertainty mapping, transfers, rendezvous, and coverage. Many of these calculations suit ground tools as well as onboard use. They are optional for a basic orbit predictor or ADCS loop.

## Select an analysis chain

| Task | Possible selection |
| --- | --- |
| Convert geometric observations | `tinyobservations` with explicitly transformed states |
| Estimate velocity from three positions | Gibbs or Herrick–Gibbs, chosen for suitable arc geometry |
| Initialize from three optical directions | Gauss, Laplace or the supported Double-r branch |
| Refine model parameters | `tinybatchls` with an application prediction callback |
| Fit conventional SGP4 mean elements | `tinysgp4fit` from a nearby GP seed and TEME observations |
| Map uncertainty | Linearized `tinycovariance` or a `tinyunscented` moment transform |
| Estimate maneuvers/access | Selected transfer, relative-motion or coverage helpers |

Initial orbit estimation and differential correction are different tasks. A Gauss/Laplace estimate may need more observations and refinement: in the synthetic comparison sample, generating-truth differences reached **2.77 km** for Gauss and **169 km** for finite-difference Laplace despite close agreement with equivalent implementations. Those are sample outcomes, not general error bounds.

Lambert supports zero-revolution short/long-way paths; Double-r has an elliptic minor-arc profile; Hohmann helpers assume circular endpoint orbits. Hill/CW assumes a circular chief and small relative separation. Edelbaum is an ideal quasi-circular low-thrust estimate, not general steering. Coverage scanning can miss short/grazing events between samples. The module links below are the contracts to use when choosing a method.

## Storage and validation

Generic batch least squares needs **7,360 bytes** of caller workspace; SGP4 fitting needs **10,000 bytes** on the tested ABI, including its solver/satellite scratch. Stack, callbacks and runtime are additional. These algorithms have explicit iteration/work limits but may be unsuitable for a fast control loop.

Vallado comparisons pin [CelesTrak/fundamentals-of-astrodynamics, `49df0479c950917c0a2c5ac31f2db6477bba57f6`](https://github.com/CelesTrak/fundamentals-of-astrodynamics/tree/49df0479c950917c0a2c5ac31f2db6477bba57f6). Additional independent equivalents include SciPy, NumPy, FilterPy, Skyfield, pymap3d, GeographicLib and python-sgp4 2.27. Module READMEs record exact versions, pinned sources, tolerances, sampled differences and documented oracle discrepancies. No external oracle is a runtime or ordinary-test dependency.

## Algorithm glossary

| Algorithm | Purpose and supported profile | Inputs → outputs; connections |
| --- | --- | --- |
| [tinyobservations](tinyobservations/README.md) | Relative RA/Dec/range/rates, fixed-station azimuth/elevation/range/rates, and coordinate helpers. | Properly framed relative states ↔ geometric observables; no light time/refraction/aberration. |
| [tinygibbs](tinygibbs/README.md) | Gibbs middle velocity estimate. | Three approximately coplanar inertial positions → velocity; small-angle geometry can be ill-conditioned. |
| [tinyherrickgibbs](tinyherrickgibbs/README.md) | Herrick–Gibbs short-arc velocity approximation. | Three positions and time differences → middle velocity; alternative to Gibbs for suitable small arcs. |
| [tinygauss](tinygauss/README.md) | Classical three-angle Gauss initial orbit estimate. | Sight lines, observer positions, times and optional radius hint → middle state/ranges; eighth-degree roots and truncated f/g. |
| [tinylaplace](tinylaplace/README.md) | Classical Laplace angles-only estimate using LOS curvature. | Direction derivatives and observer kinematics → middle state; convenience three-sample helper requires equal spacing. |
| [tinydoubler](tinydoubler/README.md) | Double-r two-radius Newton/conic-flight-time solution. | Three sight lines/observers, times and radius seeds → state/ranges; elliptic e<0.95, positive far-ray minor-arc profile. |
| [tinybatchls](tinybatchls/README.md) | Weighted bounded Levenberg–Marquardt differential correction. | Callback, measurements, sigmas and seed → fitted parameters/local covariance; 1–8 parameters, ≤64 scalars. |
| [tinysgp4fit](tinysgp4fit/README.md) | Local SGP4 mean-element fit, with optional B*. | Nearby GP seed and ≤10 TEME Cartesian epochs → fitted elements/local covariance; fixed epoch/metadata, positions or positions+velocities. |
| [tinycovariance](tinycovariance/README.md) | First-order covariance mapping, numerical Jacobians, fixed-frame/orbital-coordinate wrappers. | Covariance and transformation → covariance; dimensions 1–6, explicit singular-coordinate handling. |
| [tinyunscented](tinyunscented/README.md) | Scaled Merwe sigma points and Euclidean nonlinear moment transform. | Mean/covariance and callback → transformed moments; dimensions 1–6, not a complete UKF. |
| [tinylambert](tinylambert/README.md) | Universal-variable zero-revolution two-body boundary-value solve. | Endpoint positions and flight time → departure/arrival velocities; short/long way, no multirevolution or Battin branch. |
| [tinyhohmann](tinyhohmann/README.md) | Hohmann, bielliptic and one-tangent impulsive transfer estimates. | Circular endpoint orbits and selected geometry → burn magnitudes/coast times; no finite-thrust guidance. |
| [tinyplanechange](tinyplanechange/README.md) | Inclination/node/combined impulsive plane changes. | Orbit planes, speed/flight-path angle → intersection geometry and delta-v; unchanged speed magnitude. |
| [tinyhills](tinyhills/README.md) | CW/Hill relative propagation, transition matrix and zero-position rendezvous targeting. | Circular-chief mean motion and small rotating-frame relative state → relative state/required velocity. |
| [tinyrendezvous](tinyrendezvous/README.md) | Circular coplanar phasing, moving-target Kepler/Lambert intercept, ideal Edelbaum estimate. | Orbit/transfer inputs → maneuver/coast estimates; Edelbaum is a quasi-circular low-thrust estimate, not a general steering solver. |
| [tinycoverage](tinycoverage/README.md) | Spherical geodesics, nadir-cone footprints, ENU look angles, ideal repeat tracks and bounded access scanning. | Geometry or margin callback → coverage/access; sampled scanning can miss short or grazing events. |

## Build an individual algorithm

From the collection root:

```sh
cmake -S mission/tinylambert -B build/tinylambert -DCMAKE_BUILD_TYPE=Release
cmake --build build/tinylambert
ctest --test-dir build/tinylambert --output-on-failure
```

Keep the full checkout for shared CMake helpers and dependencies across groups. Link `tinylambert::tinylambert`; public include paths and C API names are unchanged. Building one algorithm includes its required dependencies, not every algorithm in this group. See the [collection build options](../README.md#build-and-test).
