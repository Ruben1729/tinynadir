/* SPDX-License-Identifier: MIT */
#ifndef TINYGEODESY_H
#define TINYGEODESY_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum { TINYGEODESY_OK, TINYGEODESY_INVALID, TINYGEODESY_NO_CONVERGENCE } tinygeodesy_status;
typedef struct { double semi_major_axis_m, flattening; } tinygeodesy_ellipsoid;
typedef struct { double latitude_rad, longitude_rad, height_m; } tinygeodesy_coordinate;
typedef enum { TINYGEODESY_ENU, TINYGEODESY_NED, TINYGEODESY_SEZ } tinygeodesy_frame;
#define TINYGEODESY_WGS84_A_M 6378137.0
#define TINYGEODESY_WGS84_F (1.0/298.257223563)
/* Ellipsoid: a>0, 0<=f<=0.1. Coordinates in radians/metres.
   Inversion requires |position|>=a/4, excluding the nonunique central region.
   Forward heights >=-a/2. Pole longitude from inversion is zero on the exact axis.
   All failures preserve outputs. */
tinygeodesy_status tinygeodesy_to_ecef(const tinygeodesy_ellipsoid *ellipsoid,
                                     const tinygeodesy_coordinate *geodetic,double ecef_m[3]);
tinygeodesy_status tinygeodesy_from_ecef(const tinygeodesy_ellipsoid *ellipsoid,
                                       const double ecef_m[3],tinygeodesy_coordinate *geodetic);
/* Row-major rotation maps ECEF vector coordinates to the selected local frame. */
tinygeodesy_status tinygeodesy_local_frame(double latitude_rad,double longitude_rad,
                                          tinygeodesy_frame frame,double rotation[9]);
/* Geocentric latitude corresponding to a geodetic latitude on the ellipsoid surface. */
tinygeodesy_status tinygeodesy_geocentric_latitude(double geodetic_latitude_rad,
                                                 double flattening,double *latitude_rad);
#ifdef __cplusplus
}
#endif
#endif
