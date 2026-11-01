#ifndef TINYEPHEMERIS_H
#define TINYEPHEMERIS_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct {double time_s,position_m[3],velocity_m_s[3];} tinyephemeris_sample;
typedef enum {TINYEPHEMERIS_OK=0,TINYEPHEMERIS_INVALID,TINYEPHEMERIS_RANGE,TINYEPHEMERIS_ARITHMETIC} tinyephemeris_status;
/* Cubic Hermite position interpolation; returned velocity is its derivative.
   2..4096 strictly time-sorted samples in ONE frame/origin/time scale.
   Caller owns const samples; no extrapolation; failed outputs preserved.
   position/velocity outputs must be distinct and may not overlap the table. */
tinyephemeris_status tinyephemeris_evaluate(const tinyephemeris_sample *samples,size_t n,double time_s,double position_m[3],double velocity_m_s[3]);
#ifdef __cplusplus
}
#endif
#endif
