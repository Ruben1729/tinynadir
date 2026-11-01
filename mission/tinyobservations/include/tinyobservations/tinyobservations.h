/* SPDX-License-Identifier: MIT */
#ifndef TINYOBSERVATIONS_H
#define TINYOBSERVATIONS_H
#include "tinygeodesy/tinygeodesy.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef enum { TINYOBSERVATIONS_OK,TINYOBSERVATIONS_INVALID,TINYOBSERVATIONS_SINGULAR } tinyobservations_status;
typedef struct {double position_m[3],velocity_m_s[3];} tinyobservations_state;
typedef struct {double range_m,right_ascension_rad,declination_rad,range_rate_m_s,right_ascension_rate_rad_s,declination_rate_rad_s;} tinyobservations_radec;
typedef struct {double range_m,azimuth_rad,elevation_rad,range_rate_m_s,azimuth_rate_rad_s,elevation_rate_rad_s;} tinyobservations_azel;
/* Cartesian input is relative to the desired observer and in one named frame.
   Angles radians, rates rad/s. Right ascension returned in [0,2*pi).
   A direction within 1e-14 rad of a spherical pole is singular.
   No light-time, aberration, refraction, or Earth rotation model is implicit.
   All failures preserve outputs. */
tinyobservations_status tinyobservations_to_radec(const tinyobservations_state *relative,tinyobservations_radec *output);
tinyobservations_status tinyobservations_from_radec(const tinyobservations_radec *angles,tinyobservations_state *output);
/* Target ECEF velocity is the derivative of ECEF coordinates. Observer is fixed
   at the supplied geodetic station. Azimuth is clockwise from North toward East. */
tinyobservations_status tinyobservations_to_azel(const tinygeodesy_ellipsoid *ellipsoid,const tinygeodesy_coordinate *station,const tinyobservations_state *target_ecef,tinyobservations_azel *output);
tinyobservations_status tinyobservations_from_azel(const tinygeodesy_ellipsoid *ellipsoid,const tinygeodesy_coordinate *station,const tinyobservations_azel *angles,tinyobservations_state *target_ecef);
/* Generic same-frame relative state and equatorial->ecliptic vector rotation.
   Use negative obliquity for the inverse ecliptic->equatorial rotation. */
tinyobservations_status tinyobservations_relative(const tinyobservations_state *target,const tinyobservations_state *observer,tinyobservations_state *output);
tinyobservations_status tinyobservations_ecliptic(const tinyobservations_state *equatorial,double obliquity_rad,tinyobservations_state *ecliptic);
#ifdef __cplusplus
}
#endif
#endif
