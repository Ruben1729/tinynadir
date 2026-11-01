/* SPDX-License-Identifier: MIT */
#include "tinydisturbance/tinydisturbance.h"
#include "tinylinalg/tinylinalg.h"
#include <math.h>
#include <stddef.h>
#include <string.h>
static int finite_array(const double *v,size_t n){if(!v)return 0;for(size_t i=0;i<n;++i)if(!isfinite(v[i]))return 0;return 1;}
tinydisturbance_status tinydisturbance_gravity_gradient(double mu,const double position[3],const double inertia[9],double torque[3]){
    double unit[3],momentum[3]={0},cross[3],r,coefficient,lower[9],scratch[9];
    if(!torque||!isfinite(mu)||mu<=0||!finite_array(position,3)||!finite_array(inertia,9))return TINYDISTURBANCE_INVALID_ARGUMENT;
    if(tinylinalg_cholesky(3,inertia,lower,scratch,9)!=TINYLINALG_OK)return TINYDISTURBANCE_INVALID_ARGUMENT;
    if(tinylinalg_unit3(position,unit)!=TINYLINALG_OK)return TINYDISTURBANCE_DEGENERATE;
    r=tinylinalg_norm3(position);if(!isfinite(r))return TINYDISTURBANCE_ARITHMETIC_ERROR;
    coefficient=3*((mu/r)/r)/r;
    for(size_t i=0;i<3;++i)for(size_t j=0;j<3;++j)momentum[i]+=inertia[i*3+j]*unit[j];
    if(tinylinalg_cross3(unit,momentum,cross)!=TINYLINALG_OK)return TINYDISTURBANCE_ARITHMETIC_ERROR;
    for(size_t i=0;i<3;++i){cross[i]*=coefficient;if(!isfinite(cross[i]))return TINYDISTURBANCE_ARITHMETIC_ERROR;}
    memcpy(torque,cross,sizeof cross);return TINYDISTURBANCE_OK;
}
tinydisturbance_status tinydisturbance_dipole(const double dipole[3],const double field[3],double torque[3]){
    if(!torque||!finite_array(dipole,3)||!finite_array(field,3))return TINYDISTURBANCE_INVALID_ARGUMENT;
    return tinylinalg_cross3(dipole,field,torque)==TINYLINALG_OK?TINYDISTURBANCE_OK:TINYDISTURBANCE_ARITHMETIC_ERROR;
}
tinydisturbance_status tinydisturbance_force(const double arm[3],const double force[3],double torque[3]){
    if(!torque||!finite_array(arm,3)||!finite_array(force,3))return TINYDISTURBANCE_INVALID_ARGUMENT;
    return tinylinalg_cross3(arm,force,torque)==TINYLINALG_OK?TINYDISTURBANCE_OK:TINYDISTURBANCE_ARITHMETIC_ERROR;
}
tinydisturbance_status tinydisturbance_aerodynamic(double rho,double cd,double area,const double normal[3],
    const double velocity[3],const double arm[3],double force[3],double torque[3]){
    double n[3],v[3],f[3]={0},t[3],speed,exposure;
    if(!force||!torque||!isfinite(rho)||rho<0||!isfinite(cd)||cd<0||!isfinite(area)||area<0||!finite_array(normal,3)||!finite_array(velocity,3)||!finite_array(arm,3))return TINYDISTURBANCE_INVALID_ARGUMENT;
    if(tinylinalg_unit3(normal,n)!=TINYLINALG_OK)return TINYDISTURBANCE_DEGENERATE;
    speed=tinylinalg_norm3(velocity);if(!isfinite(speed))return TINYDISTURBANCE_ARITHMETIC_ERROR;
    if(speed>0){
        for(size_t i=0;i<3;++i)v[i]=velocity[i]/speed;
        exposure=fmax(0,tinylinalg_dot3(n,v));
        if(exposure>0&&rho>0&&cd>0&&area>0)for(size_t i=0;i<3;++i){f[i]=-.5*rho*cd*area*exposure*speed*velocity[i];if(!isfinite(f[i]))return TINYDISTURBANCE_ARITHMETIC_ERROR;}
    }
    if(tinylinalg_cross3(arm,f,t)!=TINYLINALG_OK)return TINYDISTURBANCE_ARITHMETIC_ERROR;
    memcpy(force,f,sizeof f);memcpy(torque,t,sizeof t);return TINYDISTURBANCE_OK;
}
tinydisturbance_status tinydisturbance_solar_pressure(double pressure,double area,double specular,double diffuse,
    const double normal[3],const double sun[3],const double arm[3],double force[3],double torque[3]){
    double n[3],s[3],f[3],t[3],cosine,scale;
    if(!force||!torque||!isfinite(pressure)||pressure<0||!isfinite(area)||area<0||!isfinite(specular)||specular<0||!isfinite(diffuse)||diffuse<0||specular+diffuse>1||!finite_array(normal,3)||!finite_array(sun,3)||!finite_array(arm,3))return TINYDISTURBANCE_INVALID_ARGUMENT;
    if(tinylinalg_unit3(normal,n)!=TINYLINALG_OK||tinylinalg_unit3(sun,s)!=TINYLINALG_OK)return TINYDISTURBANCE_DEGENERATE;
    cosine=fmax(0,tinylinalg_dot3(n,s));scale=cosine>0&&pressure>0&&area>0?-pressure*area*cosine:0;
    for(size_t i=0;i<3;++i){f[i]=scale*((1-specular)*s[i]+2*(specular*cosine+diffuse/3)*n[i]);if(!isfinite(f[i]))return TINYDISTURBANCE_ARITHMETIC_ERROR;}
    if(tinylinalg_cross3(arm,f,t)!=TINYLINALG_OK)return TINYDISTURBANCE_ARITHMETIC_ERROR;
    memcpy(force,f,sizeof f);memcpy(torque,t,sizeof t);return TINYDISTURBANCE_OK;
}
