# tinyplanechange

Part of [Mission](../README.md) in [TinyNadir](../../README.md).

Heap-free C99 ideal impulsive inclination-only, node-only, and combined inclination/node plane changes. Include `tinyplanechange/tinyplanechange.h`, link `tinyplanechange::tinyplanechange`. Speed and delta-v are metres/second; all angles radians.

The model preserves the radial velocity component and turns the transverse component `v*cos(flight_path_angle)`. It assumes one instantaneous burn at the common line of intersection and unchanged speed magnitude. It does not combine altitude changes, thrust duration, or arbitrary burn locations. Inclinations are in [0,pi]; flight path angle is in [-pi/2,pi/2].

`tinyplanechange_combined` returns required delta-v magnitude, angle between angular-momentum directions, and initial/final argument of latitude for one intersection direction. The antipodal direction is equivalent geometrically. For parallel/antiparallel normals the line is not unique; `unique_intersection=0` and both reported arguments are zero. `tinyplanechange_node` specializes to equal inclinations. Invalid/nonfinite inputs/results preserve outputs. No allocation or mutable global state is used.

Validation: **36 combined and 27 inclination-only frozen Vallado Python cases** at commit `49df0479c950917c0a2c5ac31f2db6477bba57f6`, spanning three speeds, signed flight-path angles, small/large plane changes, and varied inclinations. Intersection arguments are compared modulo pi. Maximum sampled delta-v difference was **2.28e-12 m/s** on AppleClang 21. Zero-change/node-wrapper behavior and invalid-input preservation are also tested. Raw reference outputs/generator stay ignored under `.validation/tinyplanechange/`; normal builds/tests are offline.

Independent implementation from orbital-plane vector geometry and the transverse-velocity chord equation. No AGPL source was translated. See LICENSE/NOTICE.
