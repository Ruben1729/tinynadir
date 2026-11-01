#ifndef TINYEXPONENTIAL_H
#define TINYEXPONENTIAL_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYEXPONENTIAL_OK=0,TINYEXPONENTIAL_INVALID,TINYEXPONENTIAL_RANGE,TINYEXPONENTIAL_ARITHMETIC} tinyexponential_status;
typedef struct {double base_altitude_m,base_density_kg_m3,scale_height_m;} tinyexponential_layer;
/* n=1..128 sorted layers, positive scale heights/nonnegative densities.
   Each layer is active until the next base; maximum_altitude is inclusive.
   No extrapolation below first base or above maximum; failed out preserved. */
tinyexponential_status tinyexponential_density(const tinyexponential_layer *layers,size_t n,double maximum_altitude_m,double altitude_m,double *density_kg_m3);
#ifdef __cplusplus
}
#endif
#endif
