#ifndef TINYEOP_H
#define TINYEOP_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
#define TINYEOP_MAX_SAMPLES 1024
/* Store UT1-TAI rather than UT1-UTC to avoid interpolating across leap jumps.
   Field order is xp,yp (rad), UT1-TAI(s), LOD(s), dX,dY(rad). */
typedef struct {double mjd_utc;double value[6];} tinyeop_record;
typedef struct {double xp_rad,yp_rad,dut1_s,lod_s,dx_rad,dy_rad;} tinyeop_result;
typedef struct {double mjd_tt;double x_rad,y_rad,s_rad;} tinyeop_xys_record;
typedef enum {TINYEOP_OK=0,TINYEOP_INVALID,TINYEOP_RANGE,TINYEOP_TABLE} tinyeop_status;
typedef enum {TINYEOP_LINEAR=0,TINYEOP_CUBIC=1} tinyeop_method;
tinyeop_status tinyeop_sample(const tinyeop_record *records,size_t count,double mjd_utc,double tai_minus_utc_s,tinyeop_method method,tinyeop_result *out);
tinyeop_status tinyeop_sample_xys(const tinyeop_xys_record *records,size_t count,double mjd_tt,tinyeop_method method,double xys_rad[3]);
#ifdef __cplusplus
}
#endif
#endif
