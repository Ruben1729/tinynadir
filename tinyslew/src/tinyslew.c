/* SPDX-License-Identifier: MIT */
#include "tinyslew/tinyslew.h"
#include "tinyrotation/tinyrotation.h"
#include <math.h>
#include <stddef.h>
#include <string.h>
tinyslew_status tinyslew_initialize(tinyslew_plan *plan,const double initial[4],const double target[4],double maximum_rate,double maximum_acceleration){
    tinyslew_plan next;double inverse[4],delta[4],triangle_time;tinyrotation_status status;
    if(!plan||!isfinite(maximum_rate)||maximum_rate<=0||!isfinite(maximum_acceleration)||maximum_acceleration<=0)return TINYSLEW_INVALID_ARGUMENT;
    status=tinyrotation_normalize(initial,next.initial);if(status==TINYROTATION_DEGENERATE)return TINYSLEW_DEGENERATE;if(status!=TINYROTATION_OK)return TINYSLEW_INVALID_ARGUMENT;
    status=tinyrotation_normalize(target,next.target);if(status==TINYROTATION_DEGENERATE)return TINYSLEW_DEGENERATE;if(status!=TINYROTATION_OK)return TINYSLEW_INVALID_ARGUMENT;
    if(tinyrotation_inverse(next.initial,inverse)!=TINYROTATION_OK||tinyrotation_multiply(next.target,inverse,delta)!=TINYROTATION_OK)return TINYSLEW_ARITHMETIC_ERROR;
    if(delta[0]<0)for(size_t i=0;i<4;++i)delta[i]=-delta[i];
    if(tinyrotation_to_axis_angle(delta,next.axis,&next.angle)!=TINYROTATION_OK||tinyrotation_multiply(delta,next.initial,next.target)!=TINYROTATION_OK)return TINYSLEW_ARITHMETIC_ERROR;
    next.acceleration=maximum_acceleration;
    if(next.angle==0){next.peak_rate=next.accel_time=next.cruise_time=next.duration=0;}
    else{
        triangle_time=sqrt(next.angle)/sqrt(maximum_acceleration);
        if(triangle_time<=maximum_rate/maximum_acceleration){next.accel_time=triangle_time;next.peak_rate=maximum_acceleration*triangle_time;next.cruise_time=0;}
        else{next.accel_time=maximum_rate/maximum_acceleration;next.peak_rate=maximum_rate;next.cruise_time=next.angle/maximum_rate-next.accel_time;}
        next.duration=2*next.accel_time+next.cruise_time;
        if(!isfinite(next.duration)||next.duration<=0||!isfinite(next.peak_rate))return TINYSLEW_ARITHMETIC_ERROR;
    }
    *plan=next;return TINYSLEW_OK;
}
tinyslew_status tinyslew_sample(const tinyslew_plan *plan,double t,double q[4],double rate[3],double acceleration[3]){
    double position,speed,accel,delta[4],out[4],w[3],a[3];
    if(!plan||!q||!rate||!acceleration||!isfinite(t)||t<0||!isfinite(plan->duration)||plan->duration<0)return TINYSLEW_INVALID_ARGUMENT;
    if(t>=plan->duration){if(tinyrotation_normalize(plan->target,out)!=TINYROTATION_OK)return TINYSLEW_INVALID_ARGUMENT;speed=accel=0;}
    else{
        if(!isfinite(plan->angle)||!isfinite(plan->acceleration)||!isfinite(plan->accel_time)||!isfinite(plan->cruise_time)||!isfinite(plan->peak_rate))return TINYSLEW_INVALID_ARGUMENT;
        if(t<plan->accel_time){speed=plan->acceleration*t;position=.5*speed*t;accel=plan->acceleration;}
        else if(t<plan->accel_time+plan->cruise_time){speed=plan->peak_rate;position=.5*plan->peak_rate*plan->accel_time+speed*(t-plan->accel_time);accel=0;}
        else{double remaining=plan->duration-t;speed=plan->acceleration*remaining;position=plan->angle-.5*speed*remaining;accel=-plan->acceleration;}
        if(tinyrotation_from_axis_angle(plan->axis,position,delta)!=TINYROTATION_OK||tinyrotation_multiply(delta,plan->initial,out)!=TINYROTATION_OK)return TINYSLEW_ARITHMETIC_ERROR;
    }
    for(size_t i=0;i<3;++i){if(!isfinite(plan->axis[i]))return TINYSLEW_INVALID_ARGUMENT;w[i]=-plan->axis[i]*speed;a[i]=-plan->axis[i]*accel;if(!isfinite(w[i])||!isfinite(a[i]))return TINYSLEW_ARITHMETIC_ERROR;}
    memcpy(q,out,sizeof out);memcpy(rate,w,sizeof w);memcpy(acceleration,a,sizeof a);return TINYSLEW_OK;
}
