#ifndef TINYSPACEWEATHER_H
#define TINYSPACEWEATHER_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYSPACEWEATHER_OK=0,TINYSPACEWEATHER_INVALID,TINYSPACEWEATHER_RANGE,TINYSPACEWEATHER_ARITHMETIC} tinyspaceweather_status;
typedef struct {double time_s,f107_sfu,f107_average_sfu,ap;} tinyspaceweather_sample;
/* 2..4096 sorted samples, explicit common time scale; linear interpolation
   of each field. No extrapolation or automatic 3-hour/daily Ap aggregation. */
tinyspaceweather_status tinyspaceweather_interpolate(const tinyspaceweather_sample *samples,size_t n,double time_s,tinyspaceweather_sample *out);
/* Standard 28 Kp/Ap nodes; linear interpolation between Kp thirds. */
tinyspaceweather_status tinyspaceweather_kp_to_ap(double kp,double *ap);
tinyspaceweather_status tinyspaceweather_ap_to_kp(double ap,double *kp);
#ifdef __cplusplus
}
#endif
#endif
