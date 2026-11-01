# tinyigrf

Part of [Environment](../README.md) in [TinyNadir](../../README.md).

Independent heap-free C99 synthesis of **IGRF-14**, using the official **2025.0 main field and 2025–2030 secular variation**. The supported epoch interval is **2025.0 through 2030.0 inclusive**, with maximum degree selectable per call from 1 through 13. This small current-mission coefficient profile intentionally excludes historical epochs; an unsupported date returns `TINYIGRF_RANGE` rather than extrapolating.

Public functions in `tinyigrf/tinyigrf.h` accept geocentric spherical, Cartesian ECEF, or WGS-84 geodetic positions. Distances are metres, angles radians, and magnetic outputs **tesla**. Spherical outputs are radial/outward, colatitude/southward, and longitude/eastward. Cartesian output is ECEF; rotate it to a chosen inertial reference frame before matching it to a body-frame magnetometer measurement. It is the Earth's main magnetic-field model, not an attitude estimator or a model of spacecraft magnetic interference.

```c
#include "tinyigrf/tinyigrf.h"
double magnetic_ecef[3];
/* Geodetic latitude=45 degrees, longitude=0, altitude=400 km. */
tinyigrf_status status=tinyigrf_geodetic(2026.5,0.7853981633974483,0,400000,13,magnetic_ecef);
```

The implementation evaluates associated Legendre polynomials and their derivatives with explicit Schmidt seminormalization. A polar limiting expression avoids division by zero. Maximum degree/work and storage are bounded: two 14×14 work arrays plus small vectors/scalars, no heap or mutable globals; the 105 four-double coefficient records occupy 3360 bytes of read-only data. Degree truncation reduces computation but does not reduce this profile's compiled coefficient storage. A lower degree is an approximation whose model error must be evaluated separately.

## Verification

Independent reference: **ppigrf 2.1.0**, NumPy 2.5.4, its IGRF14 coefficient table. The comparison sampled ground/400 km/other LEO/GEO radii, longitudes, colatitudes near both poles and the equator, epochs spanning the model interval, and degrees 1, 3, 8, and 13. **218 finite spherical reference cases** and **five independent geodetic→ECEF cases** are frozen under `tests/data/`. Maximum observed spherical component difference was **2.91039×10⁻¹¹ nT**; static absolute tolerance is **0.00005 nT (5×10⁻¹⁴ T)**. The geodetic-wrapper allowance is **0.001 nT**: ppigrf uses a different finite-precision geodetic conversion and its largest sampled difference was 0.000205472 nT. These are implementation differences, not physical field accuracy guarantees.

For fractional-year comparisons, the reference datetime is chosen to produce the same linear 2025→2030 coefficient interpolation fraction; this avoids confusing pandas elapsed-calendar-time interpolation with a decimal-year input. Exact north-pole reference results are excluded because ppigrf divides by zero there; this implementation's exact-pole result is separately tested against the near-pole limit. Tests also cover model limits, aliasing, nonfinite inputs, and failure preservation.

Raw reference code/packages, independent outputs, warnings, and the source-data SHA256 manifest remain in ignored `.validation/tinyigrf/`. Normal CMake/CTest builds need only the static C fixtures.

## Sources and licensing

Equations and coefficients: [NOAA/IAGA IGRF](https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field), official [IGRF-14 coefficient text](https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt). Coefficient release: [IAGA (2024), IGRF-14, DOI 10.5281/zenodo.14012302](https://doi.org/10.5281/zenodo.14012302). Code is MIT; the coefficient data have separate **CC BY 4.0** attribution in [NOTICE](NOTICE). No ppigrf source is incorporated.
