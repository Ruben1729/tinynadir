/* SPDX-License-Identifier: MIT */
#include "tinydrag/tinydrag.h"
#include <math.h>
#include <string.h>
tinydrag_status tinydrag_acceleration(double rho,const double v[3],double cd,double area,double mass,double out[3]){
    double speed,factor,result[3];
    if(!v||!out||!isfinite(rho)||rho<0||!isfinite(cd)||cd<0||!isfinite(area)||area<0||!isfinite(mass)||mass<=0)return TINYDRAG_INVALID;
    for(size_t i=0;i<3;++i)if(!isfinite(v[i]))return TINYDRAG_INVALID;
    speed=hypot(hypot(v[0],v[1]),v[2]);factor=-.5*rho*cd*(area/mass)*speed;
    for(size_t i=0;i<3;++i){result[i]=factor*v[i];if(!isfinite(result[i]))return TINYDRAG_ARITHMETIC;}
    memcpy(out,result,sizeof result);return TINYDRAG_OK;
}
