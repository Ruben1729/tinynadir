#ifndef TINYSRP_H
#define TINYSRP_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYSRP_OK=0,TINYSRP_INVALID,TINYSRP_SINGULAR,TINYSRP_ARITHMETIC} tinysrp_status;
/* Sun position and spacecraft position share an origin and frame (metres).
   Cannonball Cr, m2, kg, illumination [0,1], pressure N/m2 at reference_m.
   Result acceleration is away from Sun. Failed output preserved. */
tinysrp_status tinysrp_acceleration(const double spacecraft_m[3],const double sun_m[3],double reflectivity,double area_m2,double mass_kg,double illumination,double reference_pressure_n_m2,double reference_distance_m,double out_m_s2[3]);
#ifdef __cplusplus
}
#endif
#endif
