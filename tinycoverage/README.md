# tinycoverage

Independent heap-free C99 geometry for ground tracking and bounded access searches. The selected profiles are:

- Great-circle distance/initial azimuth and direct destination on a **sphere** of caller-selected radius. These are geocentric spherical latitudes, not a WGS-84 ellipsoidal geodesic. Coincident points use azimuth zero; an antipodal inverse with nonunique bearing returns a geometry error.
- Circular **nadir-pointed conical** sensor footprint, limited by the spherical Earth horizon. This does not approximate an arbitrary off-nadir elliptical footprint silently.
- ECEF station/object range and local geodetic ENU azimuth/elevation. The station's location and local-basis latitude/longitude must be consistent; atmosphere/refraction is excluded. Zenith azimuth is defined as zero.
- Ideal repeating-ground-track nodal period and a two-body semimajor-axis estimate. The supplied effective Earth rotation period may account for node drift, but this function does not solve a J2 repeat-ground-track design by itself.
- Callback-based access/event intervals: the callback evaluates a scalar margin (for example elevation minus a station mask), allowing any orbit propagator, Sun/Moon source, or uploaded ephemeris to be selected independently.

All geometry uses **metres, seconds, radians**, and explicit Earth/body parameters. No orbit propagator, frame model, EOP, clocks, sensor drivers, or downloaded data are hidden dependencies. Outputs are caller-owned and preserved on failure.

## Access scanning

`tinycoverage_scan` samples a caller-specified interval and bisects sampled sign changes. Its **64-interval caller workspace** is 1024 bytes on the tested ABI. Call count is explicitly bounded (maximum accepted budget 1,000,000); each crossing has at most 64 bisections. There is no heap or shared mutable state. Bounded errors identify callback failure, storage exhaustion, work exhaustion, or insufficient floating-point time resolution. Output intervals/count/call count are committed only on complete success; scratch may change.

**An interval shorter than the sample step, or a grazing event without a sampled sign change, can be missed.** The application must choose a step consistent with its event geometry and scheduling requirements. This sampling search is not a mathematical proof of continuous visibility. Interval endpoints clipped by the requested time window are returned at that boundary.

```c
#include "tinycoverage/tinycoverage.h"
tinycoverage_footprint footprint;
/* Spherical Earth, 400 km altitude, 10-degree half-angle nadir cone. */
tinycoverage_nadir_footprint(6371200,400000,0.17453292519943295,&footprint);
```

For ground access, a callback can propagate with `tinykepler` or `tinysgp4`, transform to ECEF with `tinyframes`, call `tinycoverage_observe`, and return `elevation_rad - minimum_elevation_rad`. A separate occultation/shadow predicate can instead use `tinyshadow`. These choices belong to the application.

## Verification

**264 frozen independent geometry cases** compare spherical inverse/direct results with **GeographicLib 2.1** (zero flattening), ECEF range/azimuth/elevation with **pymap3d 3.2.0**, and footprint intersections/areas with **SciPy 1.18.1** root solving and quadrature. Static limits are 0.1 micrometre for ranges, 10⁻¹² rad for angles, 10 micrometres for footprint radius, and an absolute/relative area allowance for large footprints. Four independently captured SciPy `brentq` access intervals exercise the generic callback scanner to 10⁻⁷ seconds. Tests cover capacity/callback/work-budget errors and output preservation.

Raw reference packages, independently generated geometry and interval results, and comparison scripts remain only in ignored `.validation/tinycoverage/`. Normal C tests replay `tests/data/reference.h` offline. Implementation comparisons do not imply physical accuracy for an ellipsoidal Earth or real sensor/atmospheric geometry.

MIT for this independent implementation. Reference concepts: Vallado's spherical mission geometry and standard great-circle/ENU equations; executable reference documentation: [GeographicLib](https://geographiclib.sourceforge.io/Python/doc/), [pymap3d](https://geospace-code.github.io/pymap3d/), [SciPy brentq](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html). No source from those packages is incorporated.
