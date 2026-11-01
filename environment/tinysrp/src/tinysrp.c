/* SPDX-License-Identifier: MIT */
#include "tinysrp/tinysrp.h"
#include <math.h>
#include <string.h>
tinysrp_status tinysrp_acceleration(const double r[3],const double sun[3],double cr,double area,double mass,double light,double pressure,double reference,double out[3]){
    double d[3],distance,factor,result[3],ratio;
    if(!r||!sun||!out||!isfinite(cr)||cr<0||!isfinite(area)||area<0||!isfinite(mass)||mass<=0||!isfinite(light)||light<0||light>1||!isfinite(pressure)||pressure<0||!isfinite(reference)||reference<=0)return TINYSRP_INVALID;
    for(size_t i=0;i<3;++i){if(!isfinite(r[i])||!isfinite(sun[i]))return TINYSRP_INVALID;d[i]=r[i]-sun[i];}
    distance=hypot(hypot(d[0],d[1]),d[2]);if(distance==0)return TINYSRP_SINGULAR;
    ratio=reference/distance;factor=light*pressure*cr*(area/mass)*ratio*ratio;
    for(size_t i=0;i<3;++i){result[i]=factor*(d[i]/distance);if(!isfinite(result[i]))return TINYSRP_ARITHMETIC;}
    memcpy(out,result,sizeof result);return TINYSRP_OK;
}
