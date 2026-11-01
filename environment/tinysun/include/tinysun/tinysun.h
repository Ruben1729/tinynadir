#ifndef TINYSUN_H
#define TINYSUN_H
#ifdef __cplusplus
extern "C" {
#endif
typedef struct {double position_m[3],velocity_m_s[3];} tinysun_state;
typedef enum {TINYSUN_OK=0,TINYSUN_INVALID,TINYSUN_RANGE,TINYSUN_GEOMETRY} tinysun_status;
/* Geometric Earth-to-Sun state, BCRS/ICRS-oriented axes, TDB two-part JD.
   Selected ERFA epv00 profile: 1900-2100 (|Julian centuries from J2000|<=1).
   No light-time/aberration; this is not an apparent observed Sun direction. */
tinysun_status tinysun_geocentric(double tdb1,double tdb2,tinysun_state *out);
/* Spacecraft and Sun positions must have same origin/axes and use metres. */
tinysun_status tinysun_direction(const tinysun_state *sun,const double spacecraft_position_m[3],double unit_direction[3]);
#ifdef __cplusplus
}
#endif
#endif
