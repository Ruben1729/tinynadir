/* SPDX-License-Identifier: MIT */
#include "tinygyro/tinygyro.h"
#include "tinyrotation/tinyrotation.h"
#include <math.h>
#include <stddef.h>
#include <string.h>

static int corrected(const double gyro[3],const double bias[3],double w[3]){
    if(!gyro)return 0;
    for(size_t i=0;i<3;++i){if(!isfinite(gyro[i])||(bias&&!isfinite(bias[i])))return 0;w[i]=gyro[i]-(bias?bias[i]:0);if(!isfinite(w[i]))return 0;}
    return 1;
}
tinygyro_status tinygyro_derivative(const double q[4],const double gyro[3],const double bias[3],double qdot[4]){
    double a[4],w[3],r[4];tinyrotation_status status;
    if(!qdot||!corrected(gyro,bias,w))return TINYGYRO_INVALID_ARGUMENT;
    status=tinyrotation_normalize(q,a);
    if(status==TINYROTATION_DEGENERATE)return TINYGYRO_DEGENERATE;
    if(status!=TINYROTATION_OK)return TINYGYRO_INVALID_ARGUMENT;
    r[0]=.5*(w[0]*a[1]+w[1]*a[2]+w[2]*a[3]);
    r[1]=-.5*(w[0]*a[0]+w[1]*a[3]-w[2]*a[2]);
    r[2]=-.5*(-w[0]*a[3]+w[1]*a[0]+w[2]*a[1]);
    r[3]=-.5*(w[0]*a[2]-w[1]*a[1]+w[2]*a[0]);
    for(size_t i=0;i<4;++i)if(!isfinite(r[i]))return TINYGYRO_ARITHMETIC_ERROR;
    memcpy(qdot,r,sizeof r);return TINYGYRO_OK;
}
tinygyro_status tinygyro_step(const double q[4],const double gyro[3],const double bias[3],double dt,double out[4]){
    double a[4],w[3],delta[4];tinyrotation_status status;
    if(!out||!isfinite(dt)||dt<0||!corrected(gyro,bias,w))return TINYGYRO_INVALID_ARGUMENT;
    status=tinyrotation_normalize(q,a);
    if(status==TINYROTATION_DEGENERATE)return TINYGYRO_DEGENERATE;
    if(status!=TINYROTATION_OK)return TINYGYRO_INVALID_ARGUMENT;
    for(size_t i=0;i<3;++i){w[i]*=-dt;if(!isfinite(w[i]))return TINYGYRO_ARITHMETIC_ERROR;}
    if(tinyrotation_from_rotvec(w,delta)!=TINYROTATION_OK)return TINYGYRO_ARITHMETIC_ERROR;
    if(tinyrotation_multiply(delta,a,out)!=TINYROTATION_OK)return TINYGYRO_ARITHMETIC_ERROR;
    return TINYGYRO_OK;
}
