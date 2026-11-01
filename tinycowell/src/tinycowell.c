/* SPDX-License-Identifier: MIT */
#include "tinycowell/tinycowell.h"
#include "tinyrk4/tinyrk4.h"
#include <math.h>
#include <string.h>
static int valid(const tinycowell_model *m,const tinycowell_state *s){
    if(!m||!s||!isfinite(m->mu_m3_s2)||m->mu_m3_s2<0||(!m->mu_m3_s2&&!m->force))return 0;
    for(size_t i=0;i<3;++i)if(!isfinite(s->position_m[i])||!isfinite(s->velocity_m_s[i]))return 0;return 1;
}
static int derivative(double t,const double *y,size_t n,void *context,double *dy){
    const tinycowell_model *m=(const tinycowell_model*)context;tinycowell_state s;double a[3]={0,0,0},norm,factor=0;(void)n;
    for(size_t i=0;i<3;++i){s.position_m[i]=y[i];s.velocity_m_s[i]=y[3+i];}
    if(m->force){for(size_t i=0;i<3;++i)a[i]=NAN;if(m->force(t,&s,m->context,a))return 1;}
    if(m->mu_m3_s2){norm=hypot(hypot(y[0],y[1]),y[2]);if(!isfinite(norm)||norm==0)return 1;factor=-(m->mu_m3_s2/norm)/norm/norm;}
    for(size_t i=0;i<3;++i){dy[i]=y[3+i];dy[3+i]=factor*y[i]+a[i];if(!isfinite(dy[3+i]))return 1;}return 0;
}
static void unpack(const double *y,tinycowell_state *out){tinycowell_state result;for(size_t i=0;i<3;++i){result.position_m[i]=y[i];result.velocity_m_s[i]=y[3+i];}*out=result;}
tinycowell_status tinycowell_fixed(const tinycowell_model *m,double start,double end,const tinycowell_state *s,double max_step,size_t max_steps,tinycowell_state *out,double *w,size_t count){
    size_t steps;double span,step,*y=w;
    if(!valid(m,s)||!out||!w||count<48||!isfinite(start)||!isfinite(end)||!isfinite(max_step)||max_step<=0||!max_steps||max_steps>1000000)return TINYCOWELL_INVALID;
    span=end-start;if(!isfinite(span)||fabs(span)/max_step>(double)max_steps)return TINYCOWELL_INTEGRATION_ERROR;
    steps=(size_t)ceil(fabs(span)/max_step);if(steps==0){*out=*s;return TINYCOWELL_OK;}step=span/(double)steps;
    for(size_t i=0;i<3;++i){y[i]=s->position_m[i];y[3+i]=s->velocity_m_s[i];}
    for(size_t i=0;i<steps;++i){tinyrk4_status status=tinyrk4_step(derivative,(void*)m,6,start+(double)i*step,step,y,y,w+6,count-6);if(status!=TINYRK4_OK)return status==TINYRK4_DERIVATIVE_ERROR?TINYCOWELL_FORCE_ERROR:TINYCOWELL_INTEGRATION_ERROR;}
    unpack(y,out);return TINYCOWELL_OK;
}
tinycowell_status tinycowell_adaptive(const tinycowell_model *m,double start,double end,const tinycowell_state *s,const tinyrkf45_options *options,tinycowell_state *out,tinyrkf45_stats *stats,double *w,size_t count){
    double y[6],result[6];tinyrkf45_stats local;tinyrkf45_status status;
    if(!valid(m,s)||!out||!stats||!w||count<54)return TINYCOWELL_INVALID;
    for(size_t i=0;i<3;++i){y[i]=s->position_m[i];y[3+i]=s->velocity_m_s[i];}
    status=tinyrkf45_integrate(derivative,(void*)m,6,start,end,y,options,result,&local,w,count);
    if(status!=TINYRKF45_OK)return status==TINYRKF45_INVALID?TINYCOWELL_INVALID:(status==TINYRKF45_DERIVATIVE_ERROR?TINYCOWELL_FORCE_ERROR:TINYCOWELL_INTEGRATION_ERROR);
    unpack(result,out);*stats=local;return TINYCOWELL_OK;
}
