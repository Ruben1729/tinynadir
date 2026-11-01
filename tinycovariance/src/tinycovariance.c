/* SPDX-License-Identifier: MIT */
#include "tinycovariance/tinycovariance.h"
#include <math.h>
#include <float.h>
#define TC_PI 3.1415926535897932384626433832795
static tinycovariance_status check_psd(size_t n,const double *P)
{
    double d[6],L[36]={0},tol=128*DBL_EPSILON;size_t i,j,k;
    for(i=0;i<n;++i){if(!isfinite(P[n*i+i])||P[n*i+i]<0)return TINYCOVARIANCE_NOT_PSD;d[i]=sqrt(P[n*i+i]);for(j=0;j<n;++j)if(!isfinite(P[n*i+j]))return TINYCOVARIANCE_INVALID;}
    for(i=0;i<n;++i)for(j=0;j<=i;++j){double value,reverse;if(d[i]==0||d[j]==0){if(P[n*i+j]!=0||P[n*j+i]!=0)return TINYCOVARIANCE_NOT_PSD;value=0;}else{value=(P[n*i+j]/d[i])/d[j];reverse=(P[n*j+i]/d[j])/d[i];if(!isfinite(value)||!isfinite(reverse)||fabs(value-reverse)>tol*fmax(1,fmax(fabs(value),fabs(reverse))))return TINYCOVARIANCE_INVALID;value=.5*value+.5*reverse;if(fabs(value)>1+tol)return TINYCOVARIANCE_NOT_PSD;}for(k=0;k<j;++k)value-=L[n*i+k]*L[n*j+k];if(i==j){if(value< -tol)return TINYCOVARIANCE_NOT_PSD;L[n*i+j]=value<=0?0:sqrt(value);}else if(L[n*j+j]==0){if(fabs(value)>tol)return TINYCOVARIANCE_NOT_PSD;}else L[n*i+j]=value/L[n*j+j];}
    return TINYCOVARIANCE_OK;
}
tinycovariance_status tinycovariance_transform(size_t n,size_t m,const double *J,const double *P,double *out)
{
    double temp[36],result[36];size_t i,j,k;tinycovariance_status status;
    if(!J||!P||!out||n<1||n>6||m<1||m>6)return TINYCOVARIANCE_INVALID;
    status=check_psd(n,P);if(status!=TINYCOVARIANCE_OK)return status;
    for(i=0;i<m*n;++i)if(!isfinite(J[i]))return TINYCOVARIANCE_INVALID;
    for(i=0;i<m;++i)for(j=0;j<n;++j){double value=0;for(k=0;k<n;++k)value+=J[n*i+k]*P[n*k+j];if(!isfinite(value))return TINYCOVARIANCE_INVALID;temp[n*i+j]=value;}
    for(i=0;i<m;++i)for(j=0;j<=i;++j){double value=0;for(k=0;k<n;++k)value+=temp[n*i+k]*J[n*j+k];if(!isfinite(value))return TINYCOVARIANCE_INVALID;if(i==j&&value<0)return TINYCOVARIANCE_NOT_PSD;result[m*i+j]=value;result[m*j+i]=value;}
    for(i=0;i<m*m;++i)out[i]=result[i];
    return TINYCOVARIANCE_OK;
}
tinycovariance_status tinycovariance_jacobian(size_t n,size_t m,const double *input,const double *step,const double *period,tinycovariance_model model,void *context,double *out)
{
    double x[6],plus[6],minus[6],J[36];size_t i,j;
    if(!input||!step||!model||!out||n<1||n>6||m<1||m>6)return TINYCOVARIANCE_INVALID;
    for(i=0;i<n;++i)if(!isfinite(input[i])||!isfinite(step[i])||step[i]<=0||!isfinite(input[i]+step[i])||!isfinite(input[i]-step[i])||input[i]+step[i]==input[i]||input[i]-step[i]==input[i])return TINYCOVARIANCE_INVALID;
    if(period)for(i=0;i<m;++i)if(!isfinite(period[i])||period[i]<0)return TINYCOVARIANCE_INVALID;
    for(j=0;j<n;++j){for(i=0;i<n;++i)x[i]=input[i];for(i=0;i<m;++i){plus[i]=NAN;minus[i]=NAN;}x[j]=input[j]+step[j];if(!model(x,n,plus,m,context))return TINYCOVARIANCE_MODEL_ERROR;x[j]=input[j]-step[j];if(!model(x,n,minus,m,context))return TINYCOVARIANCE_MODEL_ERROR;for(i=0;i<m;++i){double difference;if(!isfinite(plus[i])||!isfinite(minus[i]))return TINYCOVARIANCE_MODEL_ERROR;difference=plus[i]-minus[i];if(period&&period[i]>0)difference=remainder(difference,period[i]);J[i*n+j]=difference/(2*step[j]);if(!isfinite(J[i*n+j]))return TINYCOVARIANCE_INVALID;}}
    for(i=0;i<m*n;++i)out[i]=J[i];
    return TINYCOVARIANCE_OK;
}
tinycovariance_status tinycovariance_rotate_state(const double R[9],const double P[36],double out[36])
{
    double J[36]={0};size_t i,j,k;
    if(!R||!P||!out)return TINYCOVARIANCE_INVALID;
    for(i=0;i<9;++i)if(!isfinite(R[i]))return TINYCOVARIANCE_INVALID;
    for(i=0;i<3;++i)for(j=0;j<3;++j){double v=0;for(k=0;k<3;++k)v+=R[3*i+k]*R[3*j+k];if(fabs(v-(i==j?1.:0.))>1e-10)return TINYCOVARIANCE_INVALID;J[6*i+j]=R[3*i+j];J[6*(i+3)+j+3]=R[3*i+j];}
    return tinycovariance_transform(6,6,J,P,out);
}
typedef struct {double mu;int mode;} orbit_context;
static int orbit_model(const double *x,size_t n,double *y,size_t m,void *context)
{
    orbit_context *c=context;tinyelements_state s;tinyelements_classical e;tinyelements_equinoctial q;unsigned j;
    if(n!=6||m!=6)return 0;
    if(c->mode==0||c->mode==2){for(j=0;j<3;++j){s.position_m[j]=x[j];s.velocity_m_s[j]=x[3+j];}}
    if(c->mode==0){if(tinyelements_to_classical(c->mu,&s,&e)!=TINYELEMENTS_OK||e.eccentricity<1e-8||fabs(sin(e.inclination_rad))<1e-8)return 0;y[0]=e.p_m;y[1]=e.eccentricity;y[2]=e.inclination_rad;y[3]=e.raan_rad;y[4]=e.arg_periapsis_rad;y[5]=e.true_anomaly_rad;}
    else if(c->mode==1){e.p_m=x[0];e.eccentricity=x[1];e.inclination_rad=x[2];e.raan_rad=x[3];e.arg_periapsis_rad=x[4];e.true_anomaly_rad=x[5];if(tinyelements_from_classical(c->mu,&e,&s)!=TINYELEMENTS_OK)return 0;for(j=0;j<3;++j){y[j]=s.position_m[j];y[3+j]=s.velocity_m_s[j];}}
    else if(c->mode==2){if(tinyelements_to_equinoctial(c->mu,&s,&q)!=TINYELEMENTS_OK)return 0;y[0]=q.p_m;y[1]=q.f;y[2]=q.g;y[3]=q.h;y[4]=q.k;y[5]=q.longitude_rad;}
    else {q.p_m=x[0];q.f=x[1];q.g=x[2];q.h=x[3];q.k=x[4];q.longitude_rad=x[5];if(tinyelements_from_equinoctial(c->mu,&q,&s)!=TINYELEMENTS_OK)return 0;for(j=0;j<3;++j){y[j]=s.position_m[j];y[3+j]=s.velocity_m_s[j];}}
    return 1;
}
static tinycovariance_status orbit_cov(double mu,int mode,const double x[6],const double P[36],double y[6],double out[36])
{
    double steps[6],periods[6]={0},J[36],result[36],value[6],h=cbrt(DBL_EPSILON);size_t i;orbit_context context={mu,mode};tinycovariance_status status;
    if(!isfinite(mu)||mu<=0||!P||!out||!x||!y)return TINYCOVARIANCE_INVALID;
    if(!orbit_model(x,6,value,6,&context))return TINYCOVARIANCE_SINGULAR;
    if(mode==0||mode==2){double rp=hypot(hypot(x[0],x[1]),x[2]),vp=hypot(hypot(x[3],x[4]),x[5]);for(i=0;i<3;++i){steps[i]=h*fmax(1,rp);steps[3+i]=h*fmax(1,vp);}periods[5]=2*TC_PI;if(mode==0){periods[3]=2*TC_PI;periods[4]=2*TC_PI;}}
    else {for(i=0;i<6;++i)steps[i]=h*fmax(1,fabs(x[i]));if(mode==1){if(x[1]<1e-8||fabs(sin(x[2]))<1e-8)return TINYCOVARIANCE_SINGULAR;steps[0]=fmin(steps[0],.25*x[0]);steps[1]=fmin(steps[1],.25*x[1]);steps[2]=fmin(steps[2],.25*fmin(x[2],TC_PI-x[2]));}}
    status=tinycovariance_jacobian(6,6,x,steps,periods,orbit_model,&context,J);if(status!=TINYCOVARIANCE_OK)return status;
    status=tinycovariance_transform(6,6,J,P,result);if(status!=TINYCOVARIANCE_OK)return status;
    for(i=0;i<6;++i)y[i]=value[i];
    for(i=0;i<36;++i)out[i]=result[i];
    return TINYCOVARIANCE_OK;
}
tinycovariance_status tinycovariance_cartesian_to_classical(double mu,const tinyelements_state *s,const double P[36],tinyelements_classical *out,double Pc[36])
{
    double x[6],y[6],result[36];size_t j;tinyelements_classical e;tinycovariance_status status;if(!s||!out||!Pc)return TINYCOVARIANCE_INVALID;
    for(j=0;j<3;++j){x[j]=s->position_m[j];x[3+j]=s->velocity_m_s[j];}status=orbit_cov(mu,0,x,P,y,result);if(status!=TINYCOVARIANCE_OK)return status;
    e.p_m=y[0];e.eccentricity=y[1];e.inclination_rad=y[2];e.raan_rad=y[3];e.arg_periapsis_rad=y[4];e.true_anomaly_rad=y[5];for(j=0;j<36;++j)Pc[j]=result[j];*out=e;return TINYCOVARIANCE_OK;
}
tinycovariance_status tinycovariance_classical_to_cartesian(double mu,const tinyelements_classical *e,const double P[36],tinyelements_state *out,double Pc[36])
{
    double x[6],y[6],result[36];size_t j;tinyelements_state s;tinycovariance_status status;if(!e||!out||!Pc)return TINYCOVARIANCE_INVALID;
    x[0]=e->p_m;x[1]=e->eccentricity;x[2]=e->inclination_rad;x[3]=e->raan_rad;x[4]=e->arg_periapsis_rad;x[5]=e->true_anomaly_rad;status=orbit_cov(mu,1,x,P,y,result);if(status!=TINYCOVARIANCE_OK)return status;
    for(j=0;j<3;++j){s.position_m[j]=y[j];s.velocity_m_s[j]=y[3+j];}for(j=0;j<36;++j)Pc[j]=result[j];*out=s;return TINYCOVARIANCE_OK;
}
tinycovariance_status tinycovariance_cartesian_to_equinoctial(double mu,const tinyelements_state *s,const double P[36],tinyelements_equinoctial *out,double Pc[36])
{
    double x[6],y[6],result[36];size_t j;tinyelements_equinoctial q;tinycovariance_status status;if(!s||!out||!Pc)return TINYCOVARIANCE_INVALID;
    for(j=0;j<3;++j){x[j]=s->position_m[j];x[3+j]=s->velocity_m_s[j];}status=orbit_cov(mu,2,x,P,y,result);if(status!=TINYCOVARIANCE_OK)return status;
    q.p_m=y[0];q.f=y[1];q.g=y[2];q.h=y[3];q.k=y[4];q.longitude_rad=y[5];for(j=0;j<36;++j)Pc[j]=result[j];*out=q;return TINYCOVARIANCE_OK;
}
tinycovariance_status tinycovariance_equinoctial_to_cartesian(double mu,const tinyelements_equinoctial *q,const double P[36],tinyelements_state *out,double Pc[36])
{
    double x[6],y[6],result[36];size_t j;tinyelements_state s;tinycovariance_status status;if(!q||!out||!Pc)return TINYCOVARIANCE_INVALID;
    x[0]=q->p_m;x[1]=q->f;x[2]=q->g;x[3]=q->h;x[4]=q->k;x[5]=q->longitude_rad;status=orbit_cov(mu,3,x,P,y,result);if(status!=TINYCOVARIANCE_OK)return status;
    for(j=0;j<3;++j){s.position_m[j]=y[j];s.velocity_m_s[j]=y[3+j];}for(j=0;j<36;++j)Pc[j]=result[j];*out=s;return TINYCOVARIANCE_OK;
}
