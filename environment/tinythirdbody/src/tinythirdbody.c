/* SPDX-License-Identifier: MIT */
#include "tinythirdbody/tinythirdbody.h"
#include <math.h>
#include <string.h>
tinythirdbody_status tinythirdbody_acceleration(const double r[3],const double b[3],double mu,double out[3]){
    double distance,body_radius,q=0,f,scale,result[3],d[3];
    if(!r||!b||!out||!isfinite(mu)||mu<0)return TINYTHIRDBODY_INVALID;
    for(size_t i=0;i<3;++i){if(!isfinite(r[i])||!isfinite(b[i]))return TINYTHIRDBODY_INVALID;d[i]=b[i]-r[i];if(!isfinite(d[i]))return TINYTHIRDBODY_ARITHMETIC;}
    distance=hypot(hypot(d[0],d[1]),d[2]);body_radius=hypot(hypot(b[0],b[1]),b[2]);
    if(distance==0||body_radius==0)return TINYTHIRDBODY_SINGULAR;
    if(!isfinite(distance)||!isfinite(body_radius))return TINYTHIRDBODY_ARITHMETIC;
    for(size_t i=0;i<3;++i){double u=r[i]/body_radius,v=b[i]/body_radius;q+=u*(u-2*v);}
    f=(fabs(q)<.5)?expm1(1.5*log1p(q)):pow(distance/body_radius,3)-1;
    scale=-(mu/distance)/distance/distance;
    for(size_t i=0;i<3;++i){result[i]=scale*(r[i]+f*b[i]);if(!isfinite(result[i]))return TINYTHIRDBODY_ARITHMETIC;}
    memcpy(out,result,sizeof result);return TINYTHIRDBODY_OK;
}
