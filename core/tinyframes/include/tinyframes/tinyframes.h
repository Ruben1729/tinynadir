#ifndef TINYFRAMES_H
#define TINYFRAMES_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {
 TINYFRAMES_GCRF=0,TINYFRAMES_ITRF,TINYFRAMES_TEME,TINYFRAMES_J2000,
 TINYFRAMES_MOD1976,TINYFRAMES_TOD1980,TINYFRAMES_CIRS,TINYFRAMES_TIRS,TINYFRAMES_PEF1982
} tinyframes_frame;
typedef enum {TINYFRAMES_IAU2006A=0,TINYFRAMES_IAU2000B=1} tinyframes_model;
typedef struct {
 double tt[2],ut1[2]; /* two-part Julian dates in the named scales */
 double xp_rad,yp_rad,dx_rad,dy_rad; /* IERS polar motion / CIP corrections */
 double lod_s; /* excess length of day; used for velocity rates */
} tinyframes_epoch;
typedef enum {TINYFRAMES_OK=0,TINYFRAMES_INVALID,TINYFRAMES_RANGE} tinyframes_status;
/* Matrix transforms coordinates: output = R * input at the same epoch.
   Supported dates: |Julian centuries from J2000|<=1. TEME/PEF use GMST82
   and polar motion with s'=0. Modern GCRF/CIRS paths use selected IAU model. */
tinyframes_status tinyframes_matrix(const tinyframes_epoch *epoch,tinyframes_frame from,tinyframes_frame to,tinyframes_model model,double out[3][3]);
/* SI position and coordinate velocity. Includes matrix time derivative using
   a centered 1-second interval, holding xp/yp/dX/dY/LOD fixed. Distinct output
   position/velocity arrays required; either may alias the corresponding input. */
tinyframes_status tinyframes_transform(const tinyframes_epoch *epoch,tinyframes_frame from,tinyframes_frame to,tinyframes_model model,const double position_m[3],const double velocity_m_s[3],double out_position_m[3],double out_velocity_m_s[3]);
#ifdef __cplusplus
}
#endif
#endif
