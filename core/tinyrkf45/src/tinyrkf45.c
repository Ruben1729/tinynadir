/* SPDX-License-Identifier: MIT */
#include "tinyrkf45/tinyrkf45.h"
#include <math.h>
#include <string.h>
static int finite_values(const double *v,size_t n){if(!v)return 0;for(size_t i=0;i<n;++i)if(!isfinite(v[i]))return 0;return 1;}
static int evaluate(tinyrkf45_rhs f,void *ctx,size_t n,double t,const double *y,double *k){
    if(!isfinite(t)||!finite_values(y,n))return 0;
    for(size_t i=0;i<n;++i)k[i]=NAN;
    return f(t,y,n,ctx,k)==0&&finite_values(k,n);
}
tinyrkf45_status tinyrkf45_integrate(tinyrkf45_rhs f,void *ctx,size_t n,double start,double end,const double *state,const tinyrkf45_options *o,double *out,tinyrkf45_stats *stats,double *w,size_t count){
    static const double a[6]={0,.25,.375,12.0/13,1,.5};
    static const double b[6][5]={{0},{.25},{3.0/32,9.0/32},{1932.0/2197,-7200.0/2197,7296.0/2197},{439.0/216,-8,3680.0/513,-845.0/4104},{-8.0/27,2,-3544.0/2565,1859.0/4104,-11.0/40}};
    static const double c4[6]={25.0/216,0,1408.0/2565,2197.0/4104,-.2,0};
    static const double c5[6]={16.0/135,0,6656.0/12825,28561.0/56430,-9.0/50,2.0/55};
    double *y,*temp,*candidate,t,h,direction;tinyrkf45_stats s={0,0,0,0,0};
    if(!f||!o||!out||!stats||!w||!n||n>TINYRKF45_MAX_STATE||count<9*n||!isfinite(start)||!isfinite(end)||!finite_values(state,n)||
       !isfinite(o->absolute_tolerance)||o->absolute_tolerance<=0||!isfinite(o->relative_tolerance)||o->relative_tolerance<0||
       !isfinite(o->initial_step)||o->initial_step<=0||!isfinite(o->min_step)||o->min_step<=0||!isfinite(o->max_step)||o->max_step<o->min_step||!o->max_attempts||o->max_attempts>1000000)return TINYRKF45_INVALID;
    y=w+6*n;temp=w+7*n;candidate=w+8*n;memcpy(y,state,n*sizeof *y);t=start;direction=end>=start?1:-1;
    h=direction*fmin(o->max_step,fmax(o->min_step,o->initial_step));
    for(size_t attempt=0;attempt<o->max_attempts&&t!=end;++attempt){
        double error=0,remaining=end-t,previous_h;
        if(!isfinite(remaining))return TINYRKF45_ARITHMETIC;
        if(fabs(h)>fabs(remaining))h=remaining;
        if(t+h==t)return TINYRKF45_STEP_UNDERFLOW;
        for(size_t stage=0;stage<6;++stage){
            for(size_t i=0;i<n;++i){double dy=0;for(size_t j=0;j<stage;++j)dy+=b[stage][j]*w[j*n+i];temp[i]=y[i]+h*dy;}
            if(!evaluate(f,ctx,n,t+a[stage]*h,temp,w+stage*n))return TINYRKF45_DERIVATIVE_ERROR;
            ++s.evaluations;
        }
        for(size_t i=0;i<n;++i){
            double v4=y[i],v5=y[i],scale;
            for(size_t j=0;j<6;++j){v4+=h*c4[j]*w[j*n+i];v5+=h*c5[j]*w[j*n+i];}
            if(!isfinite(v4)||!isfinite(v5))return TINYRKF45_ARITHMETIC;
            candidate[i]=v5;scale=o->absolute_tolerance+o->relative_tolerance*fmax(fabs(y[i]),fabs(v5));
            error=fmax(error,fabs(v5-v4)/scale);
        }
        if(!isfinite(error))return TINYRKF45_ARITHMETIC;
        previous_h=h;
        if(error<=1){memcpy(y,candidate,n*sizeof *y);t=(h==remaining)?end:t+h;++s.accepted;s.last_step=h;}
        else ++s.rejected;
        h*=error==0?5:fmin(5,fmax(.2,.9*pow(error,-.2)));
        if(fabs(h)>o->max_step)h=direction*o->max_step;
        if(fabs(h)<o->min_step){if(error>1&&fabs(previous_h)<=o->min_step)return TINYRKF45_STEP_UNDERFLOW;h=direction*o->min_step;}
    }
    if(t!=end)return TINYRKF45_STEP_LIMIT;
    s.final_time=t;memcpy(out,y,n*sizeof *out);*stats=s;return TINYRKF45_OK;
}
