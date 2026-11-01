/* SPDX-License-Identifier: MIT */
#include "tinylambert/tinylambert.h"
#include "tinystumpff/tinystumpff.h"
#include <math.h>
#include <float.h>
#define TL_PI 3.1415926535897932384626433832795
/* Return 0 below physical y=0 boundary, -1 on arithmetic failure. */
static int flight_time(double z,double rsum,double A,double *time,double *y)
{
    tinystumpff_values C;double ratio,K,Q,u,sh,ch;
    if(z< -40000 || z>=4*TL_PI*TL_PI || tinystumpff_evaluate(z,&C)!=TINYSTUMPFF_OK || !(C.c[2]>0))return -1;
    u=sqrt(fabs(z));
    ratio=1.4142135623730950488*(z>=0 ? cos(u/2) : cosh(u/2));
    *y=rsum-A*ratio;
    if(!isfinite(*y))return -1;
    if(*y<=0)return 0;
    if(u<.1){K=C.c[3]/(C.c[2]*sqrt(C.c[2]));Q=1-ratio*K;}
    else if(z>0){sh=sin(u/2);ch=cos(u/2);K=(u-2*sh*ch)/(2.8284271247461900976*sh*sh*sh);Q=(2*sh-u*ch)/(2*sh*sh*sh);}
    else {sh=sinh(u/2);ch=cosh(u/2);K=(2*sh*ch-u)/(2.8284271247461900976*sh*sh*sh);Q=(u*ch-2*sh)/(2*sh*sh*sh);}
    *time=sqrt(*y)*(rsum*K+A*Q);
    return isfinite(*time) ? 1 : -1;
}
tinylambert_status tinylambert_solve(double mu,const double r1[3],const double r2[3],double dt,tinylambert_path path,tinylambert_solution *out)
{
    double a,b,scale,vscale,tscale,u[3],v[3],c,cr[3],crossnorm,A,tau,lo,hi,z,t,y,f,g,gd;
    tinylambert_solution result;unsigned i,iteration;int valid;
    if(!r1 || !r2 || !out || !isfinite(mu) || mu<=0 || !isfinite(dt) || dt<=0 ||
       (path!=TINYLAMBERT_SHORT_WAY && path!=TINYLAMBERT_LONG_WAY))return TINYLAMBERT_INVALID;
    for(i=0;i<3;++i)if(!isfinite(r1[i]) || !isfinite(r2[i]))return TINYLAMBERT_INVALID;
    a=hypot(hypot(r1[0],r1[1]),r1[2]);b=hypot(hypot(r2[0],r2[1]),r2[2]);
    if(!(a>0) || !(b>0) || !isfinite(a) || !isfinite(b))return TINYLAMBERT_INVALID;
    scale=fmax(a,b);vscale=sqrt(mu/scale);tscale=scale/vscale;
    if(!(vscale>0) || !isfinite(vscale) || !isfinite(tscale))return TINYLAMBERT_RANGE;
    for(i=0;i<3;++i){u[i]=r1[i]/a;v[i]=r2[i]/b;}
    c=u[0]*v[0]+u[1]*v[1]+u[2]*v[2];c=fmax(-1,fmin(1,c));
    cr[0]=u[1]*v[2]-u[2]*v[1];cr[1]=u[2]*v[0]-u[0]*v[2];cr[2]=u[0]*v[1]-u[1]*v[0];
    crossnorm=hypot(hypot(cr[0],cr[1]),cr[2]);
    if(crossnorm<1e-10)return TINYLAMBERT_GEOMETRY;
    a/=scale;b/=scale;A=sqrt(a*b*(1+c));if(path==TINYLAMBERT_LONG_WAY)A=-A;
    tau=dt/tscale;if(!isfinite(tau) || tau<1e-6 || tau>1e8)return TINYLAMBERT_RANGE;
    lo=0;hi=4*TL_PI*TL_PI*(1-1e-8);
    if(flight_time(0,a+b,A,&t,&y)!=1)return TINYLAMBERT_RANGE;
    if(t>tau) {
        lo=-1;
        for(iteration=0;iteration<16;++iteration) {
            valid=flight_time(lo,a+b,A,&t,&y);
            if(valid<0)return TINYLAMBERT_RANGE;
            if(!valid || t<=tau)break;
            lo*=2;
        }
        if(iteration==16)return TINYLAMBERT_RANGE;
    }
    if(flight_time(hi,a+b,A,&t,&y)!=1 || t<tau)return TINYLAMBERT_RANGE;
    z=.5*(lo+hi);
    for(iteration=0;iteration<160;++iteration) {
        valid=flight_time(z,a+b,A,&t,&y);
        if(valid<0)return TINYLAMBERT_RANGE;
        if(valid && fabs(t-tau)<=8*DBL_EPSILON*fmax(1.,tau))break;
        if(!valid || t<tau)lo=z;else hi=z;
        {double next=lo+.5*(hi-lo);if(next==z)break;z=next;}
    }
    if(iteration==160 || flight_time(z,a+b,A,&t,&y)!=1 || fabs(t-tau)>1e-10*fmax(1.,tau))return TINYLAMBERT_NO_CONVERGENCE;
    f=1-y/a;g=A*sqrt(y);gd=1-y/b;
    if(g==0 || !isfinite(g))return TINYLAMBERT_GEOMETRY;
    for(i=0;i<3;++i){result.departure_m_s[i]=((r2[i]/scale)-f*(r1[i]/scale))/g*vscale;result.arrival_m_s[i]=(gd*(r2[i]/scale)-(r1[i]/scale))/g*vscale;if(!isfinite(result.departure_m_s[i]) || !isfinite(result.arrival_m_s[i]))return TINYLAMBERT_RANGE;}
    result.iterations=iteration+1;*out=result;return TINYLAMBERT_OK;
}
