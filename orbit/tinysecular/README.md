# tinysecular

Part of [Orbit](../README.md) in [TinyNadir](../../README.md).

Heap-free C99 first-order J2 secular orbit model. It holds semimajor axis, eccentricity and inclination constant while advancing the right ascension of the ascending node, argument of periapsis and mean anomaly. The gravity symmetry pole is the inertial +z axis. Use this inexpensive model for approximate long-term nodal/apsidal trends; it does not contain orbital short-period variations.

Two explicit profiles share the nodal/apsidal rates:

| Profile | Mean anomaly rate |
| --- | --- |
| `TINYSECULAR_AVERAGED_J2` | Kepler mean motion plus the first-order J2 orbit-averaged correction |
| `TINYSECULAR_VALLADO_PKEPLER` | Kepler mean motion, matching the selected no-drag `pkepler` convention (Vallado Algorithm 66) |

Let `n=sqrt(mu/a^3)`, `p=a(1-e^2)`, and `F=1.5 J2 n (R/p)^2`. The rates are `Omega_dot=-F cos(i)`, `omega_dot=F(5cos(i)^2-1)/2`, and, for the averaged profile, `M_dot=n+F sqrt(1-e^2)(3cos(i)^2-1)/2`.

Distances, gravitational parameter and time use SI; angles use radians. Inputs must be elliptic with pericenter outside the explicit central-body equatorial radius and `|J2|<=0.01`. This bound defines an input profile, not an accuracy guarantee. Secular first-order truncation accumulates error as propagation time grows; omitted J2-squared, higher harmonics and perturbations can matter. Mean elements must belong to the chosen theory; TLE/SGP4 mean elements are not interchangeable with these elements. Circular/equatorial classical angle ambiguities retain the caller's convention. Returned angles lie in `[0,2pi)`. To limit angle-reduction loss, angular inputs and accumulated angle increments are bounded to 1e12 radians. Error returns preserve the output.

No dependencies, mutable global state or workspace are needed. Build with `cmake -S orbit/tinysecular -B build/secular`, `cmake --build build/secular`, and `ctest --test-dir build/secular`; link `tinysecular::tinysecular`.

## Verification

24 frozen independent SciPy 1.18.1 quadratures average the Gauss variational equations under Cartesian J2 acceleration over one Kepler orbit. They check all three averaged rates across LEO, MEO, prograde, polar and retrograde cases. The development-host maximum difference was 2.17e-19 rad/s; C tests allow 5e-16 rad/s.

72 additional final-angle fixtures compare the no-drag `VALLADO_PKEPLER` profile with Vallado's existing Python implementation at revision `49df0479c950917c0a2c5ac31f2db6477bba57f6` over positive and negative durations up to one day. Maximum angular difference was 1.12e-13 radians; tolerance is 2e-11 radians. Tests also verify reversibility, constant a/e/i, zero-J2 behavior and invalid inputs. Reference sources/raw results live only in ignored `.validation/tinysecular/`; offline tests use frozen numbers.
