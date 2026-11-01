/* SPDX-License-Identifier: MIT */
#include "tinypd/tinypd.h"
#include "tinyrotation/tinyrotation.h"
#include <math.h>
#include <stddef.h>
#include <string.h>
tinypd_status tinypd_torque(const double current[4],const double desired[4],const double rate[3],
    const double kp[3],const double kd[3],const double ff[3],double torque[3]){
    double inverse[4],error[4],out[3],sign;tinyrotation_status status;
    if(!rate||!kp||!kd||!torque)return TINYPD_INVALID_ARGUMENT;
    for(size_t i=0;i<3;++i)if(!isfinite(rate[i])||!isfinite(kp[i])||kp[i]<0||!isfinite(kd[i])||kd[i]<0||(ff&&!isfinite(ff[i])))return TINYPD_INVALID_ARGUMENT;
    status=tinyrotation_inverse(desired,inverse);
    if(status==TINYROTATION_DEGENERATE)return TINYPD_DEGENERATE;
    if(status!=TINYROTATION_OK)return TINYPD_INVALID_ARGUMENT;
    status=tinyrotation_multiply(current,inverse,error);
    if(status==TINYROTATION_DEGENERATE)return TINYPD_DEGENERATE;
    if(status!=TINYROTATION_OK)return TINYPD_INVALID_ARGUMENT;
    sign=error[0]<0?-1:1;
    for(size_t i=0;i<3;++i){out[i]=2*kp[i]*sign*error[i+1]-kd[i]*rate[i]+(ff?ff[i]:0);if(!isfinite(out[i]))return TINYPD_ARITHMETIC_ERROR;}
    memcpy(torque,out,sizeof out);return TINYPD_OK;
}
