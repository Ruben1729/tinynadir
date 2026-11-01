# tinyobservations

Part of [Mission](../README.md) in [TinyNadir](../../README.md).

Heap-free C99 observation geometry: Cartesian relative state ↔ range/right-ascension/declination and rates; ECEF state ↔ fixed-ground-station range/azimuth/elevation and rates; same-frame state subtraction; fixed equatorial/ecliptic vector rotation. Include `tinyobservations/tinyobservations.h`, link `tinyobservations::tinyobservations` (depends on `tinygeodesy`).

Positions/ranges are metres, velocity/range rate metres/second, angles radians, angular rates radians/second. RA and azimuth are returned in [0,2*pi); declination/elevation in [-pi/2,pi/2]. Azimuth is clockwise North toward East. ECEF target velocity must be the derivative of ECEF coordinate components; the geodetic station is fixed in that frame. Earth rotation/time/frame conversion belongs to `tinyframes` or the application. No atmospheric refraction, aberration, signal light-time, or visibility decision is implicit.

A zero range is invalid; directions within 1e-14 rad of a spherical pole are singular because azimuth/RA and their rates are not uniquely defined. Inverse conversion accepts explicitly supplied pole angles but a later forward conversion can report singularity. Failure preserves output; input/output state aliases are safe. No allocation or mutable global state is used.

Validation: **3 Skyfield 1.55 RA/Dec cases and 9 pymap3d 3.2.0 az/el cases**, with rates independently computed by mpmath 1.4.1 at 80 decimal digits. Maximum sampled az/el differences on AppleClang 21 were **1.87e-9 m**, **2.23e-16 rad**, and **9.10e-13 m/s**. Static tests also check inverse geometry, ecliptic rotation round trips, poles, zero range, and failure preservation. Raw comparisons and generators are ignored in `.validation/tinyobservations/`; ordinary tests use static arrays.

Raw Vallado `rv2radec` comparisons are retained for provenance: its pinned routine returns a different RA for the sampled vectors (for example 0.80978357 versus the geometric/Skyfield 0.40489179 rad for [7000,3000,2000] km). Those values are not frozen as correct expectations. This library follows direct spherical geometry and the independent Skyfield reference.

Implementation was independently written from spherical-coordinate derivatives and local-basis geometry; no upstream source was translated. See LICENSE/NOTICE.
