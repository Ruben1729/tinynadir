#ifndef TINYDRAG_H
#define TINYDRAG_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYDRAG_OK=0,TINYDRAG_INVALID,TINYDRAG_ARITHMETIC} tinydrag_status;
/* Relative velocity is spacecraft minus atmospheric velocity in a common
   frame. rho kg/m3, area m2, mass kg; output m/s2. Failed output preserved. */
tinydrag_status tinydrag_acceleration(double density_kg_m3,const double relative_velocity_m_s[3],double drag_coefficient,double area_m2,double mass_kg,double out_m_s2[3]);
#ifdef __cplusplus
}
#endif
#endif
