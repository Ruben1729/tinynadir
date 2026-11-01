#ifndef TINYTHIRDBODY_H
#define TINYTHIRDBODY_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYTHIRDBODY_OK=0,TINYTHIRDBODY_INVALID,TINYTHIRDBODY_SINGULAR,TINYTHIRDBODY_ARITHMETIC} tinythirdbody_status;
/* SI inertial positions relative to the same origin. Differential acceleration
   includes indirect acceleration of that origin. Failed calls preserve out. */
tinythirdbody_status tinythirdbody_acceleration(const double spacecraft_m[3],const double third_body_m[3],double mu_m3_s2,double out_m_s2[3]);
#ifdef __cplusplus
}
#endif
#endif
