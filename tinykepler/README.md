# tinykepler

Heap-free C99 universal-variable two-body orbit propagation for nonradial elliptic, parabolic, and hyperbolic trajectories. Link `tinykepler::tinykepler`; its only TinyNadir dependency is `tinystumpff`.

```c
#include "tinykepler/tinykepler.h"
tinykepler_state initial = {{7000000, 0, 0}, {0, 7546.05329, 0}};
tinykepler_state next;
tinykepler_status status = tinykepler_propagate(3.986004418e14, &initial, 60, &next);
```

Position is metres, velocity metres/second, gravitational parameter metres³/second², and time seconds. Use one inertial coordinate frame throughout. The input is an osculating state; TLE mean elements cannot be substituted directly. This model omits drag, oblateness, thrust, and third-body perturbations.

The solver nondimensionalizes the universal-variable equations, brackets the monotonic time equation, and performs safeguarded Newton updates. Bracketing and solving are each bounded at 128 iterations. Elliptic elapsed times are reduced modulo their computed two-body period; precision for very long spans remains limited by input/state/time precision. The supported nondimensional time magnitude is at most 1e8 after reduction. Nearly radial states (`|rhat cross vscaled| <= 1e-12*max(1,|vscaled|)`) are rejected to avoid collision singularities. Mathematical overflow/domain failures return a status and preserve output. Input and output may alias.

No allocation or mutable global state is used. A result occupies six doubles plus bounded automatic solver workspace; libm stack/runtime cost is platform-dependent.

Validation uses **42 frozen independent SciPy 1.18.1 DOP853 integrations**, nondimensionalized with `rtol=2.5e-14`, `atol=2.5e-15`. Cases cover eccentricities 0, .2, .7, .99, 1, 1.2, and 2; signed and zero propagation intervals through six hours. Maximum component differences on AppleClang 21 were **4.23e-6 m** position and **1.40e-9 m/s** velocity. Static regression limits are 1e-4 m and 1e-7 m/s; tests also check time reversal, specific-energy conservation, and failure preservation. No numerical integrator or Python package is required by normal builds/tests. Raw oracle outputs/generator are ignored in `.validation/tinykepler/`.

Reference: universal-variable f/g propagation as described in Vallado, *Fundamentals of Astrodynamics and Applications*. Implementation was independently written from the equations. The existing comparator is [SciPy solve_ivp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html), not another copy of this universal-variable solver. See LICENSE/NOTICE.
