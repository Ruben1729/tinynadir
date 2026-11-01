/* SPDX-License-Identifier: MIT */
#include "tinysun/tinysun.h"
#include "erfa.h"
#include <math.h>
#include <stddef.h>
tinysun_status tinysun_geocentric(double a,double b,tinysun_state *out) {
    double h[2][3],v[2][3],t;tinysun_state r;int i;
    if(!out||!isfinite(a)||!isfinite(b))return TINYSUN_INVALID;
    t=((a-2451545.0)+b)/36525.0;
    if(fabs(a)>1e8||fabs(b)>1e8||fabs(t)>1.0)return TINYSUN_RANGE;
    if(tinysun_erfa_Epv00(a,b,h,v))return TINYSUN_RANGE;
    for(i=0;i<3;i++) {r.position_m[i]=-h[0][i]*149597870700.0;r.velocity_m_s[i]=-h[1][i]*(149597870700.0/86400.0);if(!isfinite(r.position_m[i])||!isfinite(r.velocity_m_s[i]))return TINYSUN_RANGE;}
    *out=r;return TINYSUN_OK;
}
tinysun_status tinysun_direction(const tinysun_state *sun,const double position[3],double out[3]) {
    double v[3],r;int i;
    if(!sun||!position||!out)return TINYSUN_INVALID;
    for(i=0;i<3;i++){if(!isfinite(sun->position_m[i])||!isfinite(position[i]))return TINYSUN_INVALID;v[i]=sun->position_m[i]-position[i];}
    r=hypot(hypot(v[0],v[1]),v[2]);if(!isfinite(r)||r<=0)return TINYSUN_GEOMETRY;
    for(i=0;i<3;i++)v[i]/=r;
    for(i=0;i<3;i++)out[i]=v[i];
    return TINYSUN_OK;
}
