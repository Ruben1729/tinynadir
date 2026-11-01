#ifndef TINYMOON_H
#define TINYMOON_H
#ifdef __cplusplus
extern "C" {
#endif
typedef struct {double position_m[3],velocity_m_s[3];} tinymoon_state;
typedef enum {TINYMOON_OK=0,TINYMOON_INVALID,TINYMOON_RANGE,TINYMOON_GEOMETRY} tinymoon_status;
/* Geometric geocentric GCRS Moon state. TT two-part JD, 1900-2100 profile.
   ERFA moon98/Meeus model; no lunar light-time longitude correction.
   Model velocity neglects small precessing-ecliptic coupling. */
tinymoon_status tinymoon_geocentric(double tt1,double tt2,tinymoon_state *out);
tinymoon_status tinymoon_direction(const tinymoon_state *moon,const double spacecraft_position_m[3],double unit_direction[3]);
#ifdef __cplusplus
}
#endif
#endif
