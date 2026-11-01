/* SPDX-License-Identifier: MIT */
#include "tinymoon/tinymoon.h"
#include "erfa.h"
#include <math.h>
#include <stddef.h>
tinymoon_status tinymoon_geocentric(double a,double b,tinymoon_state *out) {
    double pv[2][3],t;tinymoon_state r;int i;
    if(!out||!isfinite(a)||!isfinite(b))return TINYMOON_INVALID;
    t=((a-2451545.0)+b)/36525.0;
    if(fabs(a)>1e8||fabs(b)>1e8||fabs(t)>1.0)return TINYMOON_RANGE;
    tinymoon_erfa_Moon98(a,b,pv);
    for(i=0;i<3;i++){r.position_m[i]=pv[0][i]*149597870700.0;r.velocity_m_s[i]=pv[1][i]*(149597870700.0/86400.0);if(!isfinite(r.position_m[i])||!isfinite(r.velocity_m_s[i]))return TINYMOON_RANGE;}
    *out=r;return TINYMOON_OK;
}
tinymoon_status tinymoon_direction(const tinymoon_state *moon,const double position[3],double out[3]) {
    double v[3],r;int i;
    if(!moon||!position||!out)return TINYMOON_INVALID;
    for(i=0;i<3;i++){if(!isfinite(moon->position_m[i])||!isfinite(position[i]))return TINYMOON_INVALID;v[i]=moon->position_m[i]-position[i];}
    r=hypot(hypot(v[0],v[1]),v[2]);if(!isfinite(r)||r<=0)return TINYMOON_GEOMETRY;
    for(i=0;i<3;i++)v[i]/=r;
    for(i=0;i<3;i++)out[i]=v[i];
    return TINYMOON_OK;
}
