# tinymsis

Part of [Environment](../README.md) in [TinyNadir](../../README.md).

Heap-free, reentrant C99 **NRLMSISE-00** neutral atmosphere model, adapted from the original public-domain NRL/Brodowski C distribution. This is the selected full NRLMSISE-00 altitude-input formulation, **not NRLMSIS 2.x**, and not a simple atmosphere approximation.

Supply day of year, UTC seconds of day, geodetic latitude/longitude, altitude, local solar time, F10.7 daily/81-day mean flux, and daily or seven-value Ap history. Inputs use **metres, radians, seconds**, with local solar time in **hours**. Solar flux should correspond to Earth's actual distance from the Sun, as specified by the original model. UTC, longitude, and local solar time should be mutually consistent; no clock, ephemeris, or space-weather lookup occurs internally.

Outputs are He/O/N2/O2/Ar/H/N/anomalous-O number densities in **m⁻³**, total density in **kg/m³**, and local/exospheric temperatures in **kelvin**. Choose regular mass density (`effective_drag_density=0`, GTD7) or effective drag density including anomalous oxygen (`=1`, GTD7D). Optional model switches preserve the original 0/1/2 behavior; switch 9 also accepts −1 for Ap history. Metric switch 0 must remain 1 to keep the public SI contract. Passing NULL switches selects all standard variations on.

The selected interface bounds altitude to **0–1500 km**, F10.7 to **40–400**, Ap values to **0–400**, and UTC seconds to `[0,86400)`. These are finite-work API limits, not empirical accuracy guarantees. Actual UTC leap-second instants need an application adapter. The original guidance recommends nominal 150/150/4 solar/geomagnetic values below 80 km.

```c
#include "tinymsis/tinymsis.h"
tinymsis_input in={180,43200,400000,0,0,12,150,150,{4,4,4,4,4,4,4},1};
tinymsis_workspace scratch;
tinymsis_output atmosphere;
tinymsis_status status=tinymsis_evaluate(&in,0,&scratch,&atmosphere);
```

Runtime storage formerly held in file globals is now caller-owned workspace (**632 bytes** on the tested 64-bit ABI), reset on every evaluation. Separate workspaces are independent. Original spline heap allocation is replaced by a fixed ten-element temporary array; the model's actual grids have at most five nodes. All coefficients are immutable and original coefficient-clamping side effects are evaluated locally. No heap or diagnostic I/O remains; additional bounded call-stack storage is required beyond the workspace. Outputs are preserved on error; scratch is not.

## Verification

Independent executable reference: the original Brodowski NRLMSISE-00 C distribution at commit **bc9a2feba4344e74201281e563332688a4d09cc3**, compiled separately. **398 frozen cases** cover all altitude layers including 32.5/72.5/120 km boundaries, equator/poles, seasons, longitudes/local times, solar/geomagnetic drivers, both GTD7/GTD7D mass definitions, daily/Ap-history modes, and individual disabled/cross-term model switches.

All 11 returned density/temperature quantities compare against independently captured original-model outputs. The optimized comparison observed **zero difference** in every sampled output; offline C tests allow 10⁻¹⁰ relative plus an absolute zero floor appropriate to each unit. Tests additionally cover reused/interleaved caller workspaces, default switches, invalid inputs, and output preservation. Matching the original model does not establish physical forecast accuracy.

Raw source, original executable, inputs/output, reference hashes, and comparison script are only in ignored `.validation/tinymsis/`. Normal CMake/CTest builds use static `tests/data/reference.h` and have no reference-code/Python/runtime data dependency.

## Source and license

Model: Picone, Hedin, Drob and Aikin (2002), *NRLMSISE-00 empirical model of the atmosphere: Statistical comparisons and scientific issues*, JGR 107(A12), DOI [10.1029/2002JA009430](https://doi.org/10.1029/2002JA009430). C source: [Dominik Brodowski's distribution](https://git.linta.de/~brodo/nrlmsise-00.git/). The upstream public-domain statement and original US Government material attribution/disclaimers are retained in [NOTICE](NOTICE); TinyNadir wrapper/adaptation contributions use [MIT](LICENSE). This is not an NRL-endorsed product. The later 2.x models are not incorporated or claimed.
