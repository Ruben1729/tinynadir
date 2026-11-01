/* SPDX-License-Identifier: MIT */
#include "tinyexponential/tinyexponential.h"
#include <math.h>
tinyexponential_status tinyexponential_density(const tinyexponential_layer *layers,size_t n,double maximum,double h,double *out){
    size_t selected=0;double rho;
    if(!layers||!out||!n||n>128||!isfinite(h)||!isfinite(maximum))return TINYEXPONENTIAL_INVALID;
    for(size_t i=0;i<n;++i){const tinyexponential_layer *p=&layers[i];if(!isfinite(p->base_altitude_m)||!isfinite(p->base_density_kg_m3)||p->base_density_kg_m3<0||!isfinite(p->scale_height_m)||p->scale_height_m<=0||(i&&p->base_altitude_m<=layers[i-1].base_altitude_m))return TINYEXPONENTIAL_INVALID;if(h>=p->base_altitude_m)selected=i;}
    if(maximum<layers[n-1].base_altitude_m)return TINYEXPONENTIAL_INVALID;
    if(h<layers[0].base_altitude_m||h>maximum)return TINYEXPONENTIAL_RANGE;
    rho=layers[selected].base_density_kg_m3*exp(-(h-layers[selected].base_altitude_m)/layers[selected].scale_height_m);
    if(!isfinite(rho))return TINYEXPONENTIAL_ARITHMETIC;
    *out=rho;return TINYEXPONENTIAL_OK;
}
