# tinysun

Part of [Environment](../README.md) in [TinyNadir](../../README.md).

Heap-free C99 geometric Earth-to-Sun position/velocity from ERFA's `epv00`, a compact analytical VSOP2000 solution. This is a permitted-source alternative to translating Vallado's separately licensed low-precision Sun routine; the actual model is explicitly ERFA epv00, not that Vallado formula.

Input is a split **TDB Julian date**; the selected profile rejects dates more than one Julian century from J2000 (approximately 1900–2100). Output is **metres and metres/second**, geometric Earth-to-Sun state in BCRS/ICRS-oriented axes. It is not a light-time/aberration-corrected apparent Sun direction. TT may be a practical approximation to TDB only when the application's accuracy budget permits it. `tinysun_direction` subtracts a same-frame spacecraft position and normalizes the line of sight.

```c
#include "tinysun/tinysun.h"
tinysun_state sun;
double spacecraft_m[3]={7000000,0,0},direction[3];
if(tinysun_geocentric(2451545.0,0,&sun)==TINYSUN_OK)
    tinysun_direction(&sun,spacecraft_m,direction);
```

All outputs are caller-owned and preserved on failure. The reference routine has fixed work and constant coefficient tables; no heap, downloads, or mutable shared state are used. Its analytical theory has physical model error; close agreement with another implementation does not remove that error. ERFA documents a maximum heliocentric Earth position difference of about 11.2 km against DE405 over 1900–2100; use supplied precise ephemerides for tighter budgets.

## Verification and provenance

**101 independent pyerfa 2.0.1.5 `epv00` reference states** across the supported two-century interval are frozen in `tests/data/reference.h`. Max observed compiled comparison difference was zero in position and velocity for the sampled build. Offline limits are 1 mm and 10⁻⁸ m/s; tests also verify direction normalization, invalid dates/data, degenerate vectors, and output preservation. Raw scripts/output remain only in ignored `.validation/tinysun/`. No Python or ERFA installation is needed to build/test the module.

ERFA source revision [1a8044cde5b7763295d472a6443387239127c6c8](https://github.com/liberfa/erfa/tree/1a8044cde5b7763295d472a6443387239127c6c8) is adapted only by private symbol namespacing. Wrapper: MIT; embedded ERFA source: its original BSD-style license, retained in [ERFA-LICENSE.txt](ERFA-LICENSE.txt), [NOTICE](NOTICE), and the source file. Scientific basis: Moisson and Bretagnon, *Celestial Mechanics and Dynamical Astronomy* 80 (2001), 205–213.
