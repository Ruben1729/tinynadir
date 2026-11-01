# Environment

[TinyNadir](../README.md) / Environment

Predict celestial reference vectors, magnetic fields, densities, accelerations, and visibility. These models consume position/time and explicit mission data; they do not require any particular orbit propagator or attitude estimator.

## Select providers and forces

| Need | Choices and connections |
| --- | --- |
| Spacecraft position | [Orbit propagator](../orbit/README.md), GNSS, or uploaded state |
| Sun/Moon geometry | Analytical `tinysun`/`tinymoon`, or supplied samples through `tinyephemeris` |
| Magnetic reference | `tinyigrf`, transformed into the estimator's common reference frame |
| Atmospheric density | Choose `tinyexponential`, `tinymsis`, or `tinydtm`; supply required drivers/history |
| Translational forces | Selected gravity, third-body, drag and SRP terms → Cowell callback |
| Illumination/visibility | `tinyshadow` → Sun-sensor gating, SRP, power or event analysis |

IGRF is a magnetic model, not an orbit predictor or attitude estimator. Atmospheres are alternative density providers rather than successive corrections. Include each acceleration once: a gravity field with C00=1 already includes central gravity. Coordinates, epochs and origins must be consistent before combining providers. The Sun/Moon APIs supply geometric approximations, not apparent/precise ephemerides.

IGRF uses the **2025–2030 IGRF-14 profile**, degree 1–13. Atmosphere models are explicitly **NRLMSISE-00** and **DTM2000**, not later MSIS/DTM revisions. Model date/input bounds and physical errors differ from implementation-comparison precision.

## Representative storage

| Resource | Measurement |
| --- | --- |
| IGRF coefficients | 3,360 bytes read-only; two 14×14 double arrays add 3,136 bytes of stack before other locals |
| NRLMSISE-00 workspace | 632 bytes on the tested ABI, plus stack/constants |
| DTM2000 workspace | 1,104 bytes on the tested ABI, plus stack/constants |
| Cortex-M7 `-Os` code + constants | IGRF 5,200; DTM2000 12,260; MSIS 41,728; Sun 48,580 bytes, before linking/dependencies/runtime |

Lower IGRF degree reduces computation but retains the compiled coefficient table. Full selected models can cost much more flash than a controller; measure the final selection.

Validation includes independently executed pyerfa, ppigrf, original NRLMSISE-00 and Orekit DTM2000, plus numerical/geometry references recorded per module. ERFA retains BSD-style terms; DTM2000 retains Apache-2.0; NRLMSISE-00 retains public-domain provenance; IGRF coefficients retain IAGA CC BY 4.0 attribution.

## Algorithm glossary

| Algorithm | Purpose and supported profile | Inputs → outputs; connections |
| --- | --- | --- |
| [tinysun](tinysun/README.md) | ERFA epv00/VSOP2000 geometric Earth-to-Sun state, approximately 1900–2100. | TDB → BCRS/ICRS-oriented position/velocity; subtract same-frame spacecraft position for direction. Not an apparent ephemeris. |
| [tinymoon](tinymoon/README.md) | ERFA moon98/Meeus geometric lunar state, approximately 1900–2100. | TT → labeled geocentric state; optional third-body/geometry input. |
| [tinyephemeris](tinyephemeris/README.md) | Cubic Hermite interpolation of caller-supplied position/velocity samples. | Table and epoch → celestial state; alternative to analytical Sun/Moon models, with no JPL archive reader. |
| [tinyigrf](tinyigrf/README.md) | IGRF-14 main-field synthesis with the 2025–2030 coefficient profile, degree 1–13. | Supported date/location → magnetic field in Tesla; ECEF/geocentric/WGS-84 interfaces. Feeds vector estimators and magnetic-actuator simulation. |
| [tinygravity](tinygravity/README.md) | Cartesian solid-harmonic gravity, configurable evaluation degree/order 0–32. | Position, body constants and packed coefficients → potential/acceleration; central/zonal/full harmonic choices use the same kernel. |
| [tinythirdbody](tinythirdbody/README.md) | Differential point-mass third-body gravity including indirect origin acceleration. | Spacecraft/body positions and body mu → acceleration; analytical or uploaded ephemerides are interchangeable providers. |
| [tinyexponential](tinyexponential/README.md) | Layered exponential density from an explicit supplied profile. | Altitude and profile → density; inexpensive alternative to empirical atmosphere models. |
| [tinymsis](tinymsis/README.md) | Full selected NRLMSISE-00 GTD7/GTD7D atmosphere with original switches/Ap history, 0–1500 km input profile. | Location, time and activity drivers → densities/temperatures; caller workspace replaces upstream runtime globals. Not MSIS 2.x. |
| [tinydtm](tinydtm/README.md) | Full selected DTM2000 thermosphere formulation and coefficients, 120–1500 km input profile. | Location, time, flux/Kp drivers → density/temperature/composition; alternative to MSIS/exponential, not DTM2012/2020. |
| [tinyspaceweather](tinyspaceweather/README.md) | Linear interpolation of supplied flux/Ap data and standard-node Kp↔Ap conversion. | Bounded table/query → activity inputs; application supplies required averaging windows/history. |
| [tinydrag](tinydrag/README.md) | Cannonball atmospheric drag. | Density, atmosphere-relative velocity, Cd, area/mass → acceleration; density model is a separate choice. |
| [tinysrp](tinysrp/README.md) | Cannonball inverse-square solar radiation pressure. | Sun geometry, pressure, optical/area/mass parameters and illumination → acceleration; choose an eclipse provider separately. |
| [tinyshadow](tinyshadow/README.md) | Apparent angular-disc illumination, cylindrical shadow, ellipsoid segment obstruction and bracketed event bisection. | Same-frame geometry → light fraction/visibility; one bracketed crossing per event solve, no automatic event scan. |

## Build an individual algorithm

From the collection root:

```sh
cmake -S environment/tinyigrf -B build/tinyigrf -DCMAKE_BUILD_TYPE=Release
cmake --build build/tinyigrf
ctest --test-dir build/tinyigrf --output-on-failure
```

Keep the full checkout for shared CMake helpers and dependencies across groups. Link `tinyigrf::tinyigrf`; public include paths and C API names are unchanged. Building one algorithm includes its required dependencies, not every algorithm in this group. See the [collection build options](../README.md#build-and-test).
