#ifndef TINYIGRF_H
#define TINYIGRF_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYIGRF_OK=0,TINYIGRF_INVALID,TINYIGRF_RANGE} tinyigrf_status;
/* IGRF-14, epoch 2025 through 2030 inclusive; degree 1..13. Inputs use SI
   metres/radians; geocentric colatitude is north-pole=0. Outputs Tesla:
   [radial outward, colatitude/south, longitude/east]. Radius >=6000000 m. */
tinyigrf_status tinyigrf_geocentric(double decimal_year,double radius_m,double colatitude_rad,double longitude_rad,unsigned degree,double out_tesla[3]);
/* Cartesian Earth-fixed inputs/outputs. Same epoch, degree, radius limits. */
tinyigrf_status tinyigrf_ecef(double decimal_year,const double position_m[3],unsigned degree,double field_tesla[3]);
/* WGS-84 geodetic latitude/longitude (radians), ellipsoid height (metres). */
tinyigrf_status tinyigrf_geodetic(double decimal_year,double latitude_rad,double longitude_rad,double height_m,unsigned degree,double field_ecef_tesla[3]);
#ifdef __cplusplus
}
#endif
#endif
