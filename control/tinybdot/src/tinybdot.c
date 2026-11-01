/* SPDX-License-Identifier: MIT */
#include "tinybdot/tinybdot.h"
#include <math.h>
#include <stddef.h>
#include <string.h>
tinybdot_status tinybdot_step(tinybdot_state *state,const double field[3],double dt,
    double gain,double tau,const double limit[3],double dipole[3],unsigned *limited_axes){
    tinybdot_state next;double command[3],alpha;unsigned flags=0;
    if(!state||!field||!limit||!dipole||!isfinite(dt)||dt<=0||!isfinite(gain)||gain<0||!isfinite(tau)||tau<0||(state->initialized!=0&&state->initialized!=1))return TINYBDOT_INVALID_ARGUMENT;
    for(size_t i=0;i<3;++i)if(!isfinite(field[i])||!isfinite(limit[i])||limit[i]<0||(state->initialized&&(!isfinite(state->previous_field[i])||!isfinite(state->filtered_derivative[i]))))return TINYBDOT_INVALID_ARGUMENT;
    next=*state;
    if(tau==0)alpha=1;else if(dt>tau)alpha=1/(1+tau/dt);else{double x=dt/tau;alpha=x/(1+x);}
    for(size_t i=0;i<3;++i){
        if(state->initialized){double raw=(field[i]-state->previous_field[i])/dt;next.filtered_derivative[i]=state->filtered_derivative[i]+alpha*(raw-state->filtered_derivative[i]);}
        else next.filtered_derivative[i]=0;
        next.previous_field[i]=field[i];command[i]=-gain*next.filtered_derivative[i];
        if(!isfinite(next.filtered_derivative[i])||!isfinite(command[i]))return TINYBDOT_ARITHMETIC_ERROR;
        if(fabs(command[i])>limit[i])flags|=1u<<i;
        command[i]=fmax(-limit[i],fmin(limit[i],command[i]));
    }
    next.initialized=1;*state=next;memcpy(dipole,command,sizeof command);if(limited_axes)*limited_axes=flags;return TINYBDOT_OK;
}
