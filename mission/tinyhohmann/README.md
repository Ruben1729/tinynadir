# tinyhohmann

Part of [Mission](../README.md) in [TinyNadir](../../README.md).

Heap-free C99 coplanar impulsive maneuver estimates between **circular endpoint orbits**: Hohmann, bielliptic, and one-tangent coasts. Include `tinyhohmann/tinyhohmann.h`; link `tinyhohmann::tinyhohmann` (which links `tinyanomaly`). Radii are metres, mu metres³/second², time seconds, and reported burns nonnegative metres/second magnitudes. These are design estimates, not inertial burn vectors, controllers, or orbit propagators.

`tinyhohmann_circular` returns two burn magnitudes and the half-transfer-ellipse coast time. Equal radii give zero burns and a half-period coast. `tinyhohmann_bielliptic` returns three burns and two half-ellipse coast times, requiring intermediate radius >= both endpoint radii.

`tinyhohmann_one_tangent` starts at transfer periapsis and ends on the outbound arc `(0,pi]`, or starts at apoapsis and ends on the inbound arc `(pi,2*pi)`. The endpoint geometry determines eccentricity and p; it returns transfer parameters and arrival flight-path angle. Elliptic/parabolic/hyperbolic coasts are supported where physically valid. Nonphysical intersections and nonfinite inputs/results are rejected. This first profile does not cover eccentric endpoint orbits. No heap or mutable global state is used; failures preserve output.

Validation compares **5 Hohmann, 15 bielliptic, and 6 one-tangent cases** with Vallado Python at commit `49df0479c950917c0a2c5ac31f2db6477bba57f6`. Maximum sampled burn difference was **2.73e-12 m/s** on AppleClang 21. One-tangent elapsed times are additionally frozen from **SciPy 1.18.1 quadrature** of `r(nu)^2/sqrt(mu*p)`, with maximum sampled difference **1.82e-12 s**.

The raw comparison records two reference limitations rather than reproducing them: the pinned Vallado one-tangent routine returns zero time for the sampled hyperbolic coast and negative time for the sampled apoapsis-start inbound coast. Static time expectations use independent positive-time quadrature, while its valid burn/eccentricity outputs are retained. All raw outputs and generators stay ignored under `.validation/tinyhohmann/`; normal tests require no Python or oracle. SciPy is BSD licensed.

Implementation was independently written from vis-viva, conic intersection geometry, and Kepler/Barker time equations; no AGPL source was translated. See LICENSE/NOTICE.
