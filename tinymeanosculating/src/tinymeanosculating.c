/* SPDX-License-Identifier: MIT; independent numerical averaging definition. */
#include "tinymeanosculating/tinymeanosculating.h"
#include "tinyelements/tinyelements.h"
#include "tinyanomaly/tinyanomaly.h"
#include <math.h>
#define PI 3.141592653589793238462643383279502884
static double wrap(double a){double r=fmod(a,2*PI);return r<0?r+2*PI:r;}
static void pack(const tinymeanosculating_mean *q,double x[6]){x[0]=q->semimajor_axis_m;x[1]=q->f;x[2]=q->g;x[3]=q->h;x[4]=q->k;x[5]=q->mean_longitude_rad;}
static void unpack(const double x[6],tinymeanosculating_mean *q){q->semimajor_axis_m=x[0];q->f=x[1];q->g=x[2];q->h=x[3];q->k=x[4];q->mean_longitude_rad=wrap(x[5]);}
static int coordinates_valid(const tinymeanosculating_model*m,const tinymeanosculating_mean*q){
 double e=hypot(q->f,q->g),t=hypot(q->h,q->k),v[6];int i;pack(q,v);
 for(i=0;i<6;i++)if(!isfinite(v[i]))return 0;
 return q->semimajor_axis_m>0&&e<=.3&&t<=tan(85*PI/180)&&fabs(q->mean_longitude_rad)<=1e12&&q->semimajor_axis_m*(1-e)>m->equatorial_radius_m;
}
static tinymeanosculating_status setup(const tinymeanosculating_model*m,const tinymeanosculating_options*o,tinymeanosculating_workspace*w){
 if(!m||!o||!w||!isfinite(m->mu_m3_s2)||!isfinite(m->equatorial_radius_m)||!isfinite(m->j2)||!isfinite(o->max_step_s)||!isfinite(o->inverse_tolerance))return TINYMEANOSCULATING_INVALID;
 if(m->mu_m3_s2<=0||m->equatorial_radius_m<=0||fabs(m->j2)>.005||o->samples<16||o->samples>512||o->samples%2||o->max_step_s<=0||o->max_steps<1||o->max_steps>1000000||o->max_iterations<1||o->max_iterations>16||o->inverse_tolerance<1e-12||o->inverse_tolerance>1e-4)return TINYMEANOSCULATING_RANGE;
 return TINYMEANOSCULATING_OK;
}
static tinymeanosculating_status coordinates(const tinymeanosculating_model*m,const tinycowell_state*s,tinymeanosculating_mean*q){
 tinyelements_state es;tinyelements_equinoctial eq;double e,varpi,mean,anomaly;int i;
 for(i=0;i<3;i++){es.position_m[i]=s->position_m[i];es.velocity_m_s[i]=s->velocity_m_s[i];}
 if(tinyelements_to_equinoctial(m->mu_m3_s2,&es,&eq)!=TINYELEMENTS_OK)return TINYMEANOSCULATING_RANGE;
 e=hypot(eq.f,eq.g);varpi=e>0?atan2(eq.g,eq.f):0;
 if(e>=1||tinyanomaly_from_true(e,remainder(eq.longitude_rad-varpi,2*PI),&mean,&anomaly)!=TINYANOMALY_OK)return TINYMEANOSCULATING_RANGE;
 q->semimajor_axis_m=eq.p_m/(1-e*e);q->f=eq.f;q->g=eq.g;q->h=eq.h;q->k=eq.k;q->mean_longitude_rad=wrap(varpi+mean);
 return coordinates_valid(m,q)?TINYMEANOSCULATING_OK:TINYMEANOSCULATING_RANGE;
}
static tinymeanosculating_status state_from_coordinates(const tinymeanosculating_model*m,const tinymeanosculating_mean*q,tinycowell_state*s){
 tinyelements_state es;tinyelements_equinoctial eq;tinyanomaly_solution anomaly;double e=hypot(q->f,q->g),varpi=e>0?atan2(q->g,q->f):0;int i;
 if(!coordinates_valid(m,q))return TINYMEANOSCULATING_RANGE;
 if(tinyanomaly_from_mean(e,q->mean_longitude_rad-varpi,&anomaly)!=TINYANOMALY_OK)return TINYMEANOSCULATING_RANGE;
 eq.p_m=q->semimajor_axis_m*(1-e*e);eq.f=q->f;eq.g=q->g;eq.h=q->h;eq.k=q->k;eq.longitude_rad=wrap(varpi+anomaly.true_anomaly_rad);
 if(tinyelements_from_equinoctial(m->mu_m3_s2,&eq,&es)!=TINYELEMENTS_OK)return TINYMEANOSCULATING_RANGE;
 for(i=0;i<3;i++){s->position_m[i]=es.position_m[i];s->velocity_m_s[i]=es.velocity_m_s[i];}return TINYMEANOSCULATING_OK;
}
static int j2_force(double epoch,const tinycowell_state*s,void*context,double a[3]){
 const tinymeanosculating_model*m=(const tinymeanosculating_model*)context;double r=hypot(hypot(s->position_m[0],s->position_m[1]),s->position_m[2]),z,factor;(void)epoch;
 if(!isfinite(r)||r<=m->equatorial_radius_m)return 1;
 z=s->position_m[2]/r;factor=1.5*m->j2*(m->mu_m3_s2/r)/r/r*pow(m->equatorial_radius_m/r,2);
 a[0]=factor*s->position_m[0]*(5*z*z-1);a[1]=factor*s->position_m[1]*(5*z*z-1);a[2]=factor*s->position_m[2]*(5*z*z-3);return 0;
}
static tinymeanosculating_status average(const tinymeanosculating_model*m,const tinycowell_state*initial,const tinymeanosculating_options*o,tinymeanosculating_workspace*w,tinymeanosculating_mean*out,tinymeanosculating_stats*stats){
 tinymeanosculating_mean center,q;tinycowell_model model;tinycowell_state s;double sum[6],x[6],period,n,dt;int branch,k;unsigned j;tinymeanosculating_status status;
 status=coordinates(m,initial,&center);if(status!=TINYMEANOSCULATING_OK)return status;
 n=sqrt(m->mu_m3_s2/center.semimajor_axis_m)/center.semimajor_axis_m;period=2*PI/n;dt=period/o->samples;
 if(!isfinite(period)||!isfinite(n)||n<=0||dt/o->max_step_s>(double)o->max_steps)return TINYMEANOSCULATING_BUDGET;
 pack(&center,sum);sum[5]=0; /* Longitude residual measured about center, fast phase removed. */
 model.mu_m3_s2=m->mu_m3_s2;model.force=j2_force;model.context=(void*)m;
 for(branch=-1;branch<=1;branch+=2){
  s=*initial;
  for(j=1;j<=o->samples/2;j++){
   double start=branch*(double)(j-1)*dt,end=branch*(double)j*dt,required=ceil(fabs(end-start)/o->max_step_s),weight=j==o->samples/2?.5:1;
   size_t steps;
   if(!isfinite(required)||required<1||required>(double)(o->max_steps-stats->integration_steps))return TINYMEANOSCULATING_BUDGET;
   steps=(size_t)required;
   if(tinycowell_fixed(&model,start,end,&s,o->max_step_s,steps,&s,w->integration,48)!=TINYCOWELL_OK)return TINYMEANOSCULATING_INTEGRATION;
   stats->integration_steps+=steps;
   status=coordinates(m,&s,&q);if(status!=TINYMEANOSCULATING_OK)return status;
   pack(&q,x);x[5]=remainder(q.mean_longitude_rad-center.mean_longitude_rad-n*end,2*PI);
   for(k=0;k<6;k++)sum[k]+=weight*x[k];
  }
 }
 for(k=0;k<6;k++)sum[k]/=o->samples;
 sum[5]+=center.mean_longitude_rad;unpack(sum,out);stats->averaging_passes++;
 return coordinates_valid(m,out)?TINYMEANOSCULATING_OK:TINYMEANOSCULATING_RANGE;
}
tinymeanosculating_status tinymeanosculating_from_osculating(const tinymeanosculating_model*m,const tinycowell_state*s,const tinymeanosculating_options*o,tinymeanosculating_workspace*w,tinymeanosculating_mean*out,tinymeanosculating_stats*stats){
 tinymeanosculating_status status;tinymeanosculating_mean result;tinymeanosculating_stats local={0,0,0};int i;
 if(!s||!out||!stats)return TINYMEANOSCULATING_INVALID;
 status=setup(m,o,w);if(status!=TINYMEANOSCULATING_OK)return status;
 for(i=0;i<3;i++)if(!isfinite(s->position_m[i])||!isfinite(s->velocity_m_s[i]))return TINYMEANOSCULATING_INVALID;
 status=average(m,s,o,w,&result,&local);if(status!=TINYMEANOSCULATING_OK)return status;
 *out=result;*stats=local;return TINYMEANOSCULATING_OK;
}
tinymeanosculating_status tinymeanosculating_to_osculating(const tinymeanosculating_model*m,const tinymeanosculating_mean*mean,const tinymeanosculating_options*o,tinymeanosculating_workspace*w,tinycowell_state*out,tinymeanosculating_stats*stats){
 tinymeanosculating_status status;tinymeanosculating_mean candidate,computed;tinymeanosculating_stats local={0,0,0};tinycowell_state state;double target[6],x[6],actual[6],delta[6];unsigned iteration;int k;
 if(!mean||!out||!stats)return TINYMEANOSCULATING_INVALID;
 status=setup(m,o,w);if(status!=TINYMEANOSCULATING_OK)return status;
 if(!coordinates_valid(m,mean))return TINYMEANOSCULATING_RANGE;
 candidate=*mean;pack(mean,target);target[5]=wrap(target[5]);
 for(iteration=0;iteration<o->max_iterations;iteration++){
  status=state_from_coordinates(m,&candidate,&state);if(status!=TINYMEANOSCULATING_OK)return status;
  status=average(m,&state,o,w,&computed,&local);if(status!=TINYMEANOSCULATING_OK)return status;
  pack(&computed,actual);pack(&candidate,x);local.residual=0;
  for(k=0;k<6;k++){
   delta[k]=k==5?remainder(target[k]-actual[k],2*PI):target[k]-actual[k];
   local.residual=fmax(local.residual,fabs(delta[k])/(k==0?target[0]:1));
   x[k]+=delta[k];
  }
  if(local.residual<=o->inverse_tolerance){*out=state;*stats=local;return TINYMEANOSCULATING_OK;}
  unpack(x,&candidate);
 }
 return TINYMEANOSCULATING_NO_CONVERGENCE;
}
