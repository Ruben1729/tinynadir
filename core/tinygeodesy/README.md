# tinygeodesy

Part of [Core](../README.md) in [TinyNadir](../../README.md).

Heap-free C99 ellipsoidal geodetic/ECEF conversion and local ENU/NED/SEZ coordinate matrices. Include `tinygeodesy/tinygeodesy.h`, link `tinygeodesy::tinygeodesy`. A `tinygeodesy_ellipsoid` supplies semimajor axis in metres and flattening; WGS84 constants are provided. Latitude/longitude are radians and height is ellipsoidal metres, not altitude above mean sea level.

The ellipsoid profile is `a>0`, `0<=f<=0.1`. Forward conversion accepts heights at least `-a/2`. Inverse conversion excludes the nonunique central region by requiring ECEF radius `>=a/4`; its fixed-point normal iteration is bounded at 100 steps. On the exact polar axis the output longitude is zero. Pole vectors produced from finite floating-point latitude may retain a tiny horizontal component and therefore retain its longitude. Failure preserves outputs.

`tinygeodesy_local_frame` returns a row-major matrix taking ECEF vector components into ENU, NED, or SEZ. It does not include Earth rotation, geoid undulation, plate motion, or a velocity-frame derivative. `tinygeodesy_geocentric_latitude` converts surface geodetic latitude to geocentric latitude. This module is independent of orbit propagation and geomagnetic-field models.

Validation: **36 static cases**, spanning WGS84/WGS72/a sphere/an artificial f=.05 ellipsoid, poles, dateline, negative heights, LEO, and GEO. Forward positions and ENU bases were compared with pymap3d 3.2.0; inverse positions use a separate mpmath 1.4.1 90-digit solution of the ellipsoid-normal equation. Maximum differences on AppleClang 21 were **7.46e-9 m** forward, **1.12e-16 rad** inverse latitude, and **7.46e-9 m** inverse height. NED/SEZ mappings, round trips and invalid central coordinates are also tested.

The raw comparison retains pymap3d inverse results as well: its closed-form approximation differed by up to 6.34e-7 rad for sampled WGS84 GEO input and 2.50e-4 rad for the artificial high-flattening case, despite small height error. Static inverse tests therefore use the high-precision normal-root oracle rather than treating that approximation as exact. Reference generators/raw output stay ignored under `.validation/tinygeodesy/`; ordinary tests are offline C99.

Reference equations are standard geodetic ellipsoid geometry; implementation was written independently. See LICENSE/NOTICE for comparison provenance.
