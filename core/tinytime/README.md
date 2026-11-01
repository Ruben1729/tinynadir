# tinytime

Part of [Core](../README.md) in [TinyNadir](../../README.md).

Heap-free C99 calendar and split-Julian-date operations. Public declarations are in `tinytime/tinytime.h`; functions return explicit status and preserve outputs on failure. Scalar arithmetic uses `double`; no clock, file, network, leap-second database, or allocation is accessed.

`tinytime_jd` stores midnight Julian date in `first` (ending in `.5`) and a fraction in `[0,1)`. Calendar conversion uses the proleptic Gregorian calendar, years 1600–9999, and uniform 86400-second days. Calendar input `second=60` is rejected: an actual UTC leap-second instant needs an application clock adapter. UTC→TAI and UTC→UT1 require the application's applicable TAI−UTC and UT1−UTC values. TAI↔TT uses the defining 32.184-second offset. Scale labels are carried by the application/function contract rather than embedded in the two-double structure.

The TT→TDB function implements a geocentric two-term periodic approximation. It is limited to Julian centuries −1 through +1 about J2000 and is not a topocentric/relativistic clock solution. Compare the stated error against mission requirements before selecting it.

```c
#include "tinytime/tinytime.h"
tinytime_jd utc, tai, tt;
if (tinytime_calendar_to_jd(2026, 11, 1, 12, 0, 0, &utc) == TINYTIME_OK &&
    tinytime_utc_to_tai(utc, 37.0, &tai) == TINYTIME_OK) {
    tinytime_tai_to_tt(tai, &tt);
}
```

## Verification

Independent reference: **pyerfa 2.0.1.5** (ERFA/SOFA algorithms), NumPy 2.5.4. Ten frozen calendar cases compare against `erfa.cal2jd`, including Gregorian century/leap boundaries, J2000, and supported year endpoints. Scale offsets, normalization, elapsed seconds, reverse calendar conversion, invalid inputs, and output preservation are also exercised.

A 10,001-date sweep across 1900–2100 compared the TDB approximation against `erfa.dtdb` for a geocentric observer. Maximum observed approximation difference was **4.32236×10⁻⁵ s (43.224 microseconds)**. Twenty independently calculated ERFA samples are frozen under `tests/data/tdb.h` with a **0.1 ms** approximation tolerance. This is an approximation error allowance, not floating-point roundoff. Calendar fractional-day tolerance is 2×10⁻¹⁶ day.

Raw reference scripts/output and the Python reference packages reside only in ignored `.validation/tinytime/`. Normal CMake/CTest builds use checked-in static values and require no Python or reference library.

The source is an independent implementation of standard Gregorian/Julian-day arithmetic and the documented periodic time approximation. References: [IAU SOFA time-scale tools](https://www.iausofa.org/current-software) and [ERFA](https://github.com/liberfa/erfa). MIT for TinyNadir contributions; no ERFA source is incorporated in this module.
