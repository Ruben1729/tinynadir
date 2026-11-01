/* SPDX-License-Identifier: MIT */
#include "tinypid/tinypid.h"
#include "tinypd/tinypd.h"
#include "tinyrotation/tinyrotation.h"
#include <math.h>
#include <stddef.h>
#include <string.h>
tinypid_status tinypid_step(tinypid_state *state,const double current[4],const double desired[4],
    const double rate[3],const double kp[3],const double kd[3],const double ki[3],const double ff[3],
    double dt,const double imax[3],const double tmax[3],double torque[3],unsigned *limited_axes){
    double base[3],inverse[4],error[4],out[3];tinypid_state next;unsigned flags=0;tinypd_status status;
    if(!state||!ki||!imax||!tmax||!torque||!isfinite(dt)||dt<0)return TINYPID_INVALID_ARGUMENT;
    for(size_t i=0;i<3;++i)if(!isfinite(state->integral[i])||!isfinite(ki[i])||ki[i]<0||!isfinite(imax[i])||imax[i]<0||!isfinite(tmax[i])||tmax[i]<0)return TINYPID_INVALID_ARGUMENT;
    status=tinypd_torque(current,desired,rate,kp,kd,ff,base);
    if(status==TINYPD_DEGENERATE)return TINYPID_DEGENERATE;
    if(status==TINYPD_ARITHMETIC_ERROR)return TINYPID_ARITHMETIC_ERROR;
    if(status!=TINYPD_OK)return TINYPID_INVALID_ARGUMENT;
    if(tinyrotation_inverse(desired,inverse)!=TINYROTATION_OK||tinyrotation_multiply(current,inverse,error)!=TINYROTATION_OK)return TINYPID_DEGENERATE;
    next=*state;
    for(size_t i=0;i<3;++i){
        double old=fmax(-imax[i],fmin(imax[i],state->integral[i]));
        double delta=2*(error[0]<0?-1:1)*error[i+1]*dt;
        double candidate=fmax(-imax[i],fmin(imax[i],old+delta)),raw;
        if(!isfinite(delta)||!isfinite(ki[i]*candidate))return TINYPID_ARITHMETIC_ERROR;
        raw=base[i]+ki[i]*candidate;
        if(!isfinite(raw))return TINYPID_ARITHMETIC_ERROR;
        if(fabs(raw)>tmax[i]&&ki[i]>0&&((candidate>old&&raw>0)||(candidate<old&&raw<0))){candidate=old;raw=base[i]+ki[i]*candidate;}
        next.integral[i]=candidate;
        if(fabs(raw)>tmax[i])flags|=1u<<i;
        out[i]=fmax(-tmax[i],fmin(tmax[i],raw));
    }
    *state=next;memcpy(torque,out,sizeof out);if(limited_axes)*limited_axes=flags;return TINYPID_OK;
}
