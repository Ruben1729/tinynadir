# tinymoon

Part of [Environment](../README.md) in [TinyNadir](../../README.md).

Heap-free C99 geometric geocentric Moon state from ERFA's `moon98`, implementing the Meeus lunar theory and its analytical velocity. This is a permitted-source alternative to a direct port of Vallado's separately licensed lunar formula; it does not claim numerical identity with that other approximation.

Supply a split **TT Julian date**; the selected profile is limited to one Julian century on either side of J2000 (approximately 1900–2100). Position/velocity are **metres and metres/second in GCRS axes**. The model omits the lunar light-time correction and a small precessing-ecliptic velocity coupling, as described by ERFA. `tinymoon_direction` subtracts a spacecraft position expressed in the same frame and units and normalizes the line of sight.

```c
#include "tinymoon/tinymoon.h"
tinymoon_state moon;
double spacecraft_m[3]={7000000,0,0},direction[3];
if(tinymoon_geocentric(2451545.0,0,&moon)==TINYMOON_OK)
    tinymoon_direction(&moon,spacecraft_m,direction);
```

The coefficient tables/counts are immutable, work is fixed, and outputs are caller-owned and preserved on failure. No heap or shared mutable state is used. This remains an approximate lunar model: ERFA reports position differences up to about 31.7 km against ELP/MPP02 in its 1950–2100 comparison. Select precise external ephemerides when the model error exceeds the mission budget.

## Verification and provenance

**101 independently executed pyerfa 2.0.1.5 `moon98` states** spanning the selected date interval are frozen under `tests/data/`. The optimized comparison observed zero difference in sampled position and velocity; an unoptimized test build showed at most 1.13687×10⁻¹³ m/s velocity difference. Static tolerances are 1 mm and 10⁻⁸ m/s. Additional tests cover normalized spacecraft-relative directions, invalid dates/data, degenerate vectors, and failure preservation. Raw oracle code/output remain only in ignored `.validation/tinymoon/`; ordinary tests are offline C.

Source basis: Meeus, *Astronomical Algorithms*, second edition (1998), lunar theory; [ERFA revision 1a8044cde5b7763295d472a6443387239127c6c8](https://github.com/liberfa/erfa/tree/1a8044cde5b7763295d472a6443387239127c6c8). Modifications are private symbol names, read-only coefficient declarations, and a bounded integer absolute expression replacing `abs`. The wrapper is MIT; embedded ERFA code retains its BSD-style copyright/license in [ERFA-LICENSE.txt](ERFA-LICENSE.txt), [NOTICE](NOTICE), and source files.
