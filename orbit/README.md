# Orbit

[TinyNadir](../README.md) / Orbit

Predict spacecraft position and velocity, or decode the GP elements used by SGP4. Choose the propagator appropriate to the supplied initial conditions and required physics. This group does not require attitude estimation or control.

## Choose a propagation path

| Starting point | A suitable path |
| --- | --- |
| Conventional TLE, including Alpha-5 | `tinytle` → `tinysgp4` |
| SGP4-compatible OMM KVN or decoded fields | `tinyomm` → SGP4 initialization/propagation |
| Osculating Cartesian state; ideal two-body coast | `tinykepler` |
| Compatible elliptic mean elements; approximate J2 trends | `tinypkepler` or `tinysecular` |
| Cartesian state with explicitly selected forces | `tinycowell` + [RK4/RKF45](../core/README.md) + chosen [environment models](../environment/README.md) |

**TLE mean elements are fitted for SGP4.** Passing them to Kepler or Cowell as osculating elements changes their meaning. The numerical J2 mean conversion has its own centered-one-orbit definition and is not a Brouwer, DSST, or SGP4 conversion. Conventional SGP4 already includes its deep-space branch; SGP4-XP is outside this collection.

For ground tracking, convert SGP4 TEME km/km/s to SI, then use [time/frames/geodesy](../core/README.md) and optional [observation helpers](../mission/README.md). GNSS/uploaded states can replace propagation in an application.

## Resource and validation notes

| Resource | Measurement, excluding additional runtime/stack |
| --- | --- |
| SGP4 satellite | 1,560 bytes on a 32-bit ABI; 1,568 on the tested 64-bit ABI |
| SGP4 code/read-only data | 17,024 bytes, Arm GCC 16.2 `-O2` measurement |
| Numerical J2 mean conversion workspace | 384 bytes; integrating one or more orbits can be expensive |

The [SGP4 provenance table](tinysgp4/README.md#compared-implementations-and-data-provenance) records pinned official Vallado C++ and aholinch revisions and **20,930 captured propagation requests**, including expected failures. The [fixed-point investigation](tinysgp4/README.md#fixed-point-investigation) records coefficient loss/error growth in the evaluated Q19.44 design; maintained code uses binary64. Other modules document their independent numerical comparisons and supported model profiles individually.

## Algorithm glossary

| Algorithm | Purpose and supported profile | Inputs → outputs; connections |
| --- | --- | --- |
| [tinytle](tinytle/README.md) | Bounded traditional/Alpha-5 TLE parsing, checksums, and decoded-element SGP4 initialization. | Two lines or validated GP elements → elements/satellite; supported ephemeris type 0. |
| [tinyomm](tinyomm/README.md) | SGP4-compatible OMM KVN parser and decoded-field adapter. | Bounded KVN/fields → GP elements/satellite; Earth, TEME, UTC, SGP4 metadata checked. JSON/XML/CSV decoding belongs to the application. |
| [tinysgp4](tinysgp4/README.md) | Vallado-compatible conventional SGP4 with near-Earth and deep-space branches. | TLE/initialized GP elements and minutes from epoch → TEME position/velocity in km and km/s. Deep-space behavior is included; SGP4-XP is outside this collection. |
| [tinykepler](tinykepler/README.md) | Universal-variable two-body propagation for nonradial elliptic/parabolic/hyperbolic states. | Osculating Cartesian state, mu, elapsed seconds → propagated state; no perturbations. |
| [tinypkepler](tinypkepler/README.md) | Selected Vallado first-order J2 secular-precession profile. | Compatible elliptic mean elements → state; fixed p/e/i, Kepler mean-anomaly rate, no J4 or drag. |
| [tinysecular](tinysecular/README.md) | First-order averaged J2 angular rates and element advancement. | Mean elements and body parameters → rates/elements; selects averaged mean-anomaly correction or the matching `tinypkepler` convention. |
| [tinymeanosculating](tinymeanosculating/README.md) | Central+J2 centered-one-orbit numerical average and bounded local inverse. | Cartesian orbit ↔ explicitly defined mean equinoctial coordinates; e≤0.3, i≤170°. This is its own filter definition, not Brouwer, DSST, or SGP4 mean elements. |
| [tinycowell](tinycowell/README.md) | Cartesian inertial orbit integration using RK4 or RKF45. | State, explicit central mu, optional force callback → state; compose gravity, drag, SRP and third bodies as needed. |

## Build an individual algorithm

From the collection root:

```sh
cmake -S orbit/tinysgp4 -B build/tinysgp4 -DCMAKE_BUILD_TYPE=Release
cmake --build build/tinysgp4
ctest --test-dir build/tinysgp4 --output-on-failure
```

Keep the full checkout for shared CMake helpers and dependencies across groups. Link `tinysgp4::tinysgp4`; public include paths and C API names are unchanged. Building one algorithm includes its required dependencies, not every algorithm in this group. See the [collection build options](../README.md#build-and-test).
