# tinydtm

Heap-free C99 **DTM2000** thermosphere model, adapted from Orekit's Apache-2.0 scalar implementation and coefficient set. This is the full selected DTM2000 formulation, **not DTM2012 or DTM2020**, and not an exponential-atmosphere approximation. Newer DTM distributions have different licensing; they are not incorporated here.

`tinydtm_evaluate` takes day of year, altitude in metres, geodetic longitude/latitude in radians, local solar time in radians (0 h = 0 rad), instantaneous/mean F10.7, and 3-hour/24-hour Kp. The selected profile accepts altitude **120–1500 km**, F10.7 **40–400**, and Kp **0–9**; these are interface bounds, not a promise of physical accuracy at every permitted combination. Inconsistent/extreme drivers can return an arithmetic failure.

Outputs are total density in **kg/m³**, local/exospheric temperatures in **kelvin**, and H/He/O/N2/O2/N partial number/mass densities in **m⁻³** and **kg/m³**. Atmosphere-relative velocity and drag acceleration belong to other modules. Space weather and local solar time are explicit inputs; nothing is downloaded or inferred from an internal clock.

```c
#include "tinydtm/tinydtm.h"
tinydtm_input in={180,400000,0,0,0,150,150,3,3};
tinydtm_workspace scratch;
tinydtm_output atmosphere;
tinydtm_status status=tinydtm_evaluate(&in,&scratch,&atmosphere);
```

The caller workspace is reusable, needs no initial values, and is **1104 bytes** on the tested 64-bit ABI. It contains latitude/local-time intermediates and a single 97-double derivative scratch reused across constituents. Coefficients are immutable; no heap, mutable globals, file I/O, or Java runtime is required. Failure preserves the output; the workspace is scratch and may change.

## Verification

Independent original executable oracle: **Orekit 13.1.5** binary Java library with Hipparchus 4.0.3, calling its unmodified `DTM2000.getDensity`. **256 cases** sample day/season boundaries, pole/equatorial latitudes, longitudes, local solar times, 120–1500 km altitudes, and varying solar/geomagnetic activity. A separately compiled instrumentation of scalar Orekit source revision **734c814858c114824f50292e838be057c614906f** exposes its already-computed local/exospheric temperatures and species concentrations; its density is additionally checked against the unmodified binary oracle.

Nine outputs per case are frozen in `tests/data/reference.h`. Maximum observed relative difference was **1.43881×10⁻¹⁴**, and maximum temperature difference **9.09495×10⁻¹³ K**. Static tolerance is 2×10⁻¹¹ relative plus a tiny absolute zero floor. Additional tests check partial-mass summation, workspace independence/reuse, invalid inputs, and output preservation. These are implementation-comparison errors, not the empirical model's physical prediction uncertainty.

Original source/jars, reference instrumentation, raw oracle output, and source/coefficient hashes stay in ignored `.validation/tinydtm/`; ordinary builds/tests use offline C fixtures only.

## Source and license

Model reference: Bruinsma, Thuillier and Barlier (2003), *The DTM-2000 empirical thermosphere model with new data assimilation and constraints at lower boundary: accuracy and properties*, Journal of Atmospheric and Solar-Terrestrial Physics 65, 1053–1070. Implementation/data source: [Orekit](https://github.com/CS-SI/Orekit/tree/734c814858c114824f50292e838be057c614906f), **Apache License 2.0**, retained in [LICENSE](LICENSE) with relevant attribution and change notices in [NOTICE](NOTICE). This module is not represented as an MIT-only port.
