# tinyframes

Part of [Core](../README.md) in [TinyNadir](../../README.md).

Heap-free C99 Earth-orientation transforms with a privately namespaced subset of **ERFA** (BSD-style licensed, derived with permission from SOFA). The public wrapper supports GCRF/GCRS-oriented axes, ITRF, SGP4 TEME, mean J2000, mean-of-date IAU-1976, true-of-date IAU-1980, CIRS, TIRS, and legacy PEF. It uses IAU-2006/2000A or the shorter IAU-2000B model. This is an orientation/state-coordinate library, not a relativistic barycentric-to-geocentric transformation or a propagator.

Provide explicitly labeled two-part **TT and UT1 Julian dates**, IERS polar motion `xp,yp`, optional celestial-pole offsets `dX,dY` (radians), and excess length of day (seconds). The profile is bounded to one Julian century on either side of J2000. It does not obtain EOP or time offsets automatically. Use `tinyeop` or application-supplied values; zero corrections are an explicit reduced-data model, not a precision guarantee.

`tinyframes_matrix` returns `output = R * input` for coordinates at the same epoch. TEME/PEF follow the conventional GMST82 path with polar motion and s′=0. GCRF↔ITRF uses the selected modern model; mean-J2000 is distinguished from GCRF by the frame-bias rotation. MOD/TOD use their named historical precession/nutation theories. CIRS/TIRS are tied to the chosen modern model. FK4/FK5 catalog transformations are not included.

`tinyframes_transform` includes the matrix time derivative in coordinate velocities using a centered **one-second difference** (±0.5 s). TT advances uniformly; UT1 advance includes the supplied LOD. The EOP corrections are held fixed across this one-second interval. This includes rotation/precession/nutation rate effects of the chosen models, but not temporal derivatives of supplied pole corrections. Inputs/outputs are metres and metres/second. Errors preserve outputs; in-place corresponding input/output arrays are supported.

```c
#include "tinyframes/tinyframes.h"
tinyframes_epoch epoch={ {2461346.5,.5}, {2461346.5,.5-69.1/86400.0}, 0,0,0,0,0 };
double teme_r[3]={7000000,0,0},teme_v[3]={0,7500,0};
double ecef_r[3],ecef_v[3];
tinyframes_status status=tinyframes_transform(&epoch,TINYFRAMES_TEME,TINYFRAMES_ITRF,
    TINYFRAMES_IAU2006A,teme_r,teme_v,ecef_r,ecef_v);
```

The SGP4 library returns km and km/s; convert those to SI before passing its state here. Frames and time scales must describe the same physical epoch.

## Verification

Independent executable reference: **pyerfa 2.0.1.5**, NumPy 2.5.4, using ERFA transforms and independent matrix composition. **648 frozen matrices** cover all 9×9 frame pairs, two models, four dates spanning 1901–2099, and nonzero polar/CIP corrections. **56 frozen SI state cases** additionally verify the velocity derivative path. Matrix component tolerance is 5×10⁻¹⁴; state tolerances are 1 micrometre and 1 micrometre/second. Observed comparison maxima were 2.22045×10⁻¹⁶ matrix components, 9.31323×10⁻¹⁰ m, and 2.34195×10⁻¹⁰ m/s. Offline tests also check orthogonality, inverse state round trips, bad dates, nonfinite data, enum errors, and output preservation.

These are differences against the same documented astronomical theories, not measurement accuracy claims. Reference code/packages and raw comparison output remain in ignored `.validation/tinyframes/`; static fixtures are under `tests/data/` and need no runtime reference dependencies.

## Sources and license

ERFA source revision **1a8044cde5b7763295d472a6443387239127c6c8** from [liberfa/erfa](https://github.com/liberfa/erfa/tree/1a8044cde5b7763295d472a6443387239127c6c8). The embedded source subset is namespaced; its original copyright and BSD-style terms remain in each file and [ERFA-LICENSE.txt](ERFA-LICENSE.txt). The wrapper is MIT. See [NOTICE](NOTICE). The astronomical coefficient tables are read-only; the kernel uses caller outputs and bounded stack arrays, with no heap or shared mutable state. Compile/link only this module when required; its full IAU-2000A series is a deliberate flash-size/accuracy choice rather than a universally minimal approximation.
