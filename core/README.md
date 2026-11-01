# Core

[TinyNadir](../README.md) / Core

Shared mathematics, time scales, coordinate systems, and numerical integration. These pieces support other groups and can also be used independently. A dependency on one helper does not include the entire group.

## Using the foundations

- `tinyrotation` supplies attitude conventions; `tinylinalg` supplies bounded matrix operations. Neither owns an estimator or controller.
- `tinytime` and `tinyeop` supply explicit time/Earth-orientation inputs for `tinyframes`. The application provides valid leap offsets and EOP windows.
- `tinyelements`, `tinyanomaly`, and `tinystumpff` support [orbit](../orbit/README.md) and [mission](../mission/README.md) algorithms. Element conversion does not convert between incompatible mean-element theories.
- `tinyrk4` and `tinyrkf45` integrate caller derivatives: select fixed steps or bounded adaptive work. They can serve Cartesian orbit or rigid-body dynamics independently.

Calendar handling uses uniform 86400-second days and rejects `second=60`; actual UTC leap-second instants need an application adapter. Frames use explicitly labeled TT/UT1 and corrections. Modern frame models are bounded to approximately 1900–2100; zero EOP corrections are an explicit approximation.

## Memory and numerical references

All kernels require binary64 doubles. Caller workspace, bounded stack, and platform math/runtime contribute to RAM use. The full frame model includes substantial read-only series: one Cortex-M7 GCC 16.2 `-Os` compile/archive measurement gave **124,608 bytes** of frame code plus constants before linking, excluding dependencies/runtime. Select individual algorithms when this fidelity is unnecessary.

Foundational numerical comparisons use NumPy 2.5.4, SciPy 1.18.1, mpmath 1.4.1, pyerfa 2.0.1.5 and module-specific equivalents. See individual READMEs for pinned references, input bounds, tolerances and failures. Embedded ERFA code in `tinyframes` retains its BSD-style license and attribution.

## Algorithm glossary

| Algorithm | Purpose and supported profile | Inputs → outputs; connections |
| --- | --- | --- |
| [tinylinalg](tinylinalg/README.md) | Small vector/matrix operations, LU/linear solves, Cholesky, symmetric eigenpairs; dimensions 1–16. | Arrays and scratch → factors/solutions; shared estimator/allocation mathematics. |
| [tinyrotation](tinyrotation/README.md) | Quaternion, DCM, axis-angle, rotation-vector, ZYX Euler conversions, composition, and SLERP. | Rotation → another representation or rotated vector; establishes attitude conventions. |
| [tinytime](tinytime/README.md) | Gregorian calendar, split Julian dates, UTC/TAI/TT/UT1, approximate geocentric TT→TDB. | Epoch and supplied time offsets → labeled time scales; calendar years 1600–9999, uniform 86400-second days, no `second=60` input. |
| [tinyeop](tinyeop/README.md) | Bounded caller-table interpolation of Earth orientation and optional IAU X/Y/s. | Epoch, table, applicable TAI−UTC → EOP; UT1−TAI interpolation avoids leap discontinuities. |
| [tinyframes](tinyframes/README.md) | Nine Earth-orientation frames, including TEME, GCRF, ITRF, CIRS/TIRS and legacy MOD/TOD/PEF; ERFA IAU-2006/2000A or 2000B. | TT, UT1, EOP and vector/state → transformed coordinates; velocity includes matrix-rate effects. Approximately 1900–2100. |
| [tinygeodesy](tinygeodesy/README.md) | Ellipsoidal geodetic/ECEF conversion and local ENU/NED/SEZ bases. | Position and explicit ellipsoid → location/local axes; ellipsoidal height, not geoid height. |
| [tinyelements](tinyelements/README.md) | Cartesian/classical/prograde modified-equinoctial conversion; RSW, NTW and perifocal bases. | Osculating state and mu ↔ elements; classical coordinates use semilatus rectum and true anomaly. |
| [tinyanomaly](tinyanomaly/README.md) | Elliptic/hyperbolic Kepler and parabolic Barker equations, anomaly conversion. | Eccentricity/anomaly → anomaly and solve status; mathematical helper, not a state propagator. |
| [tinystumpff](tinystumpff/README.md) | Stable universal-variable Stumpff functions C0–C3. | Scalar argument → four functions; supports Kepler/Lambert. |
| [tinyinterpolation](tinyinterpolation/README.md) | Linear, bounded Lagrange, and cubic Hermite interpolation. | Samples/query → value and supported derivatives; no extrapolation. |
| [tinypolynomial](tinypolynomial/README.md) | Real-coefficient polynomial roots, degrees 1–8, including complex roots. | Coefficients and caller workspace → roots/status; bounded Aberth iteration, not certified root isolation. |
| [tinyelliptic](tinyelliptic/README.md) | Complete/incomplete real Legendre F/E, inverse E, and signed ellipse arc length. | Parameter m=k² in [0,1] and angles → integrals/arcs; distinct from Kepler equations. |
| [tinyrk4](tinyrk4/README.md) | Classical fixed-step fourth-order Runge–Kutta, 1–32 components. | State, derivative callback and one step → advanced state; physics supplied separately. |
| [tinyrkf45](tinyrkf45/README.md) | Adaptive Fehlberg 4(5), 1–32 components, explicit tolerances and attempt/step bounds. | State and derivative callback → integrated state/statistics; alternative to fixed RK4 for nonstiff problems. |

## Build an individual algorithm

From the collection root:

```sh
cmake -S core/tinyrotation -B build/tinyrotation -DCMAKE_BUILD_TYPE=Release
cmake --build build/tinyrotation
ctest --test-dir build/tinyrotation --output-on-failure
```

Keep the full checkout for shared CMake helpers and dependencies across groups. Link `tinyrotation::tinyrotation`; public include paths and C API names are unchanged. Building one algorithm includes its required dependencies, not every algorithm in this group. See the [collection build options](../README.md#build-and-test).
