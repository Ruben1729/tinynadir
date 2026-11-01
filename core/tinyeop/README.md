# tinyeop

Part of [Core](../README.md) in [TinyNadir](../../README.md).

Heap-free C99 interpolation of application-supplied Earth-orientation and optional IAU X/Y/s sample windows. `tinyeop_sample` interpolates polar motion, UT1 offset, length of day, and celestial pole corrections; `tinyeop_sample_xys` interpolates precomputed X/Y/s. APIs use modified Julian dates, seconds, and radians. Tables are immutable and caller-owned; no downloads, archive readers, or implicit Earth-orientation data are included.

Choose two-point linear or local four-point cubic **Lagrange** interpolation (not a global spline). Records must be strictly increasing and finite; count is limited to 1024 and queries must lie inside the supplied interval. Each query validates the table and has bounded work. There is no extrapolation. Outputs are preserved on failure.

To avoid interpolating a leap-second jump, `tinyeop_record.value[2]` stores **UT1−TAI**, not UT1−UTC. Supply the applicable TAI−UTC offset for the query epoch; the returned `dut1_s` is UT1−UTC. UTC leap-second timestamp encoding and obtaining valid IERS samples remain application responsibilities. Ingesting a source's DUT1 requires subtracting that sample's TAI−UTC first.

```c
#include "tinyeop/tinyeop.h"
/* Columns: xp(rad), yp(rad), UT1-TAI(s), LOD(s), dX(rad), dY(rad). */
const tinyeop_record window[]={
    {61000,{1e-6,2e-6,-37.1,.001,0,0}},
    {61001,{1.1e-6,2.1e-6,-37.101,.0011,0,0}}
};
tinyeop_result eop;
tinyeop_status status=tinyeop_sample(window,2,61000.5,37,TINYEOP_LINEAR,&eop);
```

## Verification

Independent equivalent: **NumPy 2.5.4** `interp` (linear) and `polynomial.Polynomial.fit` (four-point cubic solved in another basis). Fifty-eight frozen comparisons cover all knots, interval boundaries/interiors, six EOP channels, and the X/Y/s adapter. The six-channel fixture allowance is 2×10⁻¹³ in native units (rad or s), with tighter angular-only adapter checks. Tests cover duplicate/nonfinite data, insufficient windows, out-of-range queries, leap-offset reconstruction, and failure preservation. Synthetic smoothly varying fields make the interpolation algorithm independent of a particular IERS data release; this does not validate or supply actual Earth-orientation predictions.

Raw calculations and tables are in ignored `.validation/tinyeop/`; offline tests use `tests/data/reference.h`. Reference algorithms and data concepts: [IERS Earth orientation](https://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html) and [NumPy interpolation](https://numpy.org/doc/stable/reference/generated/numpy.interp.html). Independent implementation, MIT license; no NumPy or Vallado code is incorporated.
