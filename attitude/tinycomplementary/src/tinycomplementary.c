/* SPDX-License-Identifier: MIT */
#include "tinycomplementary/tinycomplementary.h"
#include "tinylinalg/tinylinalg.h"
#include "tinyrotation/tinyrotation.h"
#include "tinygyro/tinygyro.h"
#include <math.h>
#include <string.h>
tinycomplementary_status tinycomplementary_initialize(tinycomplementary_state *state,const double q[4],const double bias[3]){
    tinycomplementary_state next;tinyrotation_status status;
    if(!state||!bias)return TINYCOMPLEMENTARY_INVALID_ARGUMENT;
    for(size_t i=0;i<3;++i)if(!isfinite(bias[i]))return TINYCOMPLEMENTARY_INVALID_ARGUMENT;
    status=tinyrotation_normalize(q,next.q);if(status==TINYROTATION_DEGENERATE)return TINYCOMPLEMENTARY_DEGENERATE;if(status!=TINYROTATION_OK)return TINYCOMPLEMENTARY_INVALID_ARGUMENT;
    memcpy(next.bias,bias,sizeof next.bias);*state=next;return TINYCOMPLEMENTARY_OK;
}
tinycomplementary_status tinycomplementary_step(tinycomplementary_state *state,const tinycomplementary_parameters *p,
    const double gyro[3],double dt,size_t count,const double *reference,const double *body,const double *weights){
    tinycomplementary_state next;double error[3]={0},max_weight=0,total=0,corrected[3],dq[4];
    if(!state||!p||!gyro||!isfinite(dt)||dt<=0||count>TINYCOMPLEMENTARY_MAX_OBSERVATIONS||(count&&(!reference||!body))||!isfinite(p->proportional_gain)||p->proportional_gain<0||!isfinite(p->integral_gain)||p->integral_gain<0||!isfinite(p->bias_limit)||p->bias_limit<0)return TINYCOMPLEMENTARY_INVALID_ARGUMENT;
    if(tinycomplementary_initialize(&next,state->q,state->bias)!=TINYCOMPLEMENTARY_OK)return TINYCOMPLEMENTARY_INVALID_ARGUMENT;
    for(size_t i=0;i<3;++i)if(!isfinite(gyro[i]))return TINYCOMPLEMENTARY_INVALID_ARGUMENT;
    for(size_t i=0;i<count;++i){double w=weights?weights[i]:1;if(!isfinite(w)||w<0)return TINYCOMPLEMENTARY_INVALID_ARGUMENT;max_weight=fmax(max_weight,w);for(size_t j=0;j<3;++j)if(!isfinite(reference[i*3+j])||!isfinite(body[i*3+j]))return TINYCOMPLEMENTARY_INVALID_ARGUMENT;}
    if(max_weight>0){
        for(size_t i=0;i<count;++i)total+=(weights?weights[i]:1)/max_weight;
        for(size_t i=0;i<count;++i){
            double weight=(weights?weights[i]:1)/max_weight/total,r[3],b[3],predicted[3],cross[3];
            if(weight==0)continue;
            if(tinylinalg_unit3(reference+3*i,r)!=TINYLINALG_OK||tinylinalg_unit3(body+3*i,b)!=TINYLINALG_OK)return TINYCOMPLEMENTARY_DEGENERATE;
            if(tinyrotation_apply(next.q,r,predicted)!=TINYROTATION_OK||tinylinalg_cross3(b,predicted,cross)!=TINYLINALG_OK)return TINYCOMPLEMENTARY_ARITHMETIC_ERROR;
            for(size_t j=0;j<3;++j)error[j]+=weight*cross[j];
        }
    }
    for(size_t i=0;i<3;++i){
        double bias=next.bias[i]-p->integral_gain*error[i]*dt;
        if(!isfinite(bias))return TINYCOMPLEMENTARY_ARITHMETIC_ERROR;
        next.bias[i]=fmax(-p->bias_limit,fmin(p->bias_limit,bias));
        corrected[i]=gyro[i]-next.bias[i]+p->proportional_gain*error[i];
        if(!isfinite(corrected[i]))return TINYCOMPLEMENTARY_ARITHMETIC_ERROR;
    }
    if(tinygyro_derivative(next.q,corrected,NULL,dq)!=TINYGYRO_OK)return TINYCOMPLEMENTARY_ARITHMETIC_ERROR;
    for(size_t i=0;i<4;++i){dq[i]=next.q[i]+dt*dq[i];if(!isfinite(dq[i]))return TINYCOMPLEMENTARY_ARITHMETIC_ERROR;}
    if(tinyrotation_normalize(dq,next.q)!=TINYROTATION_OK)return TINYCOMPLEMENTARY_ARITHMETIC_ERROR;
    *state=next;return TINYCOMPLEMENTARY_OK;
}
