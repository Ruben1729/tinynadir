# tinyelements

Part of [Core](../README.md) in [TinyNadir](../../README.md).

Heap-free C99 conversion among Cartesian states, classical orbital elements, and prograde modified equinoctial elements. It also constructs RSW/RTN, NTW, and PQW local coordinate bases. Include `tinyelements/tinyelements.h`; link `tinyelements::tinyelements`.

All positions and semilatus recta are metres, velocities metres/second, angles radians, and `mu` metres³/second². Classical elements use **semilatus rectum p**, not semimajor axis, so elliptic/parabolic/hyperbolic cases use one representation. Forward conversion rejects an unphysical `1+e*cos(nu)<=0`. Inversion rejects zero-position/radial singularities. Successful returned angular coordinates are in `[0,2*pi)`.

For `e<1e-12`, eccentricity is reported as zero, argument of periapsis is zero, and true anomaly means argument of latitude. For equatorial states `abs(sin(i))<1e-12`, ascending node is zero and argument of periapsis is the signed longitude of periapsis. Circular equatorial states have both node/periapsis zero and signed true longitude in the anomaly field. Retrograde equatorial classical states are supported. Modified equinoctial elements use `f=e*cos(Omega+omega)`, `g=e*sin(Omega+omega)`, `h=tan(i/2)*cos(Omega)`, `k=tan(i/2)*sin(Omega)`, `L=Omega+omega+nu`; inclinations within 1e-8 rad of pi are rejected for this prograde convention.

Row-major local matrices map inertial vector components to the named local basis. Rotating a velocity vector is not the same as constructing a rotating-frame relative velocity; the latter also needs frame angular-rate terms. These conversions do not convert SGP4 mean elements into osculating classical elements. Failures preserve outputs; no heap or mutable global state is used.

Validation: **15 frozen Vallado Python cases** at commit `49df0479c950917c0a2c5ac31f2db6477bba57f6`, with eccentricities .1/.6/.99/1/1.4 and inclinations .1/1.3/2.9 rad, cover `coe2rv`, `rv2coe`, RSW/NTW basis matrices and PQW components. Maximum sampled differences on AppleClang 21 were **9.32e-10 m**, **3.64e-12 m/s**, and **1.67e-15 rad**. Modified equinoctial round trips, six circular/equatorial/retrograde cases, and failure preservation are additional invariants. Static test fixtures require no external library, downloads, or Python. Raw comparisons/generator are ignored under `.validation/tinyelements/`.

Implementation was independently written from vector/conic equations; no Vallado AGPL source was translated. The Vallado code was executed only as an oracle. See NOTICE and LICENSE.
