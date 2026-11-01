/* SPDX-License-Identifier: MIT */
#include "tinyelliptic/tinyelliptic.h"
#include <float.h>
#include <math.h>
#include <stddef.h>

#define TE_PI 3.14159265358979323846264338327950288
#define TE_HALF_PI 1.57079632679489661923132169163975144
/* Duplication followed by symmetric Taylor corrections. The argument ranges
   used below are x,y in [0,1], z=1; no arbitrary-scale RF/RD API is exposed. */
static int symmetric_rf(double x,double y,double z,double *answer){
    for(size_t iteration=0;iteration<64;++iteration){
        double mean=(x+y+z)/3,dx=(mean-x)/mean,dy=(mean-y)/mean,dz=(mean-z)/mean;
        if(fmax(fabs(dx),fmax(fabs(dy),fabs(dz)))<.0025){
            double e2=dx*dy-dz*dz,e3=dx*dy*dz;
            double series=1-e2/10+e3/14+e2*e2/24-3*e2*e3/44-5*e2*e2*e2/208+3*e3*e3/104+e2*e2*e3/16;
            *answer=series/sqrt(mean);return isfinite(*answer);
        }
        double sx=sqrt(x),sy=sqrt(y),sz=sqrt(z),lambda=sx*sy+sx*sz+sy*sz;
        x=.25*(x+lambda);y=.25*(y+lambda);z=.25*(z+lambda);
    }
    return 0;
}
static int symmetric_rd(double x,double y,double z,double *answer){
    double sum=0,factor=1;
    for(size_t iteration=0;iteration<64;++iteration){
        double mean=(x+y+3*z)/5,dx=(mean-x)/mean,dy=(mean-y)/mean,dz=(mean-z)/mean;
        if(fmax(fabs(dx),fmax(fabs(dy),fabs(dz)))<.0015){
            double ea=dx*dy,eb=dz*dz,ec=ea-eb,ed=ea-6*eb,ee=ed+2*ec;
            double series=1+ed*(-3.0/14+9.0/88*ed-9.0/52*dz*ee)+dz*(ee/6+dz*(-9.0/22*ec+dz*3.0/26*ea));
            *answer=3*sum+factor*series/(mean*sqrt(mean));return isfinite(*answer);
        }
        double sx=sqrt(x),sy=sqrt(y),sz=sqrt(z),lambda=sx*sy+sx*sz+sy*sz;
        sum+=factor/(sz*(z+lambda));factor*=.25;
        x=.25*(x+lambda);y=.25*(y+lambda);z=.25*(z+lambda);
    }
    return 0;
}
static tinyelliptic_status principal(double phi,double m,int second,double *value){
    double s=sin(phi),c=fabs(phi)==TE_HALF_PI?0:cos(phi),rf,rd=0;
    if(m==0){*value=phi;return TINYELLIPTIC_OK;}
    if(m==1){
        if(second){*value=s;return TINYELLIPTIC_OK;}
        if(fabs(phi)>=TE_HALF_PI)return TINYELLIPTIC_SINGULAR;
        *value=asinh(tan(phi));return TINYELLIPTIC_OK;
    }
    if(!symmetric_rf(c*c,1-m*s*s,1,&rf))return TINYELLIPTIC_NO_CONVERGENCE;
    if(second&&!symmetric_rd(c*c,1-m*s*s,1,&rd))return TINYELLIPTIC_NO_CONVERGENCE;
    *value=s*rf-(second?m*s*s*s*rd/3:0);
    return isfinite(*value)?TINYELLIPTIC_OK:TINYELLIPTIC_ARITHMETIC_ERROR;
}
static tinyelliptic_status extended(double phi,double m,int second,double *value){
    double n,remainder,complete,tail;tinyelliptic_status status;
    if(m==0){*value=phi;return TINYELLIPTIC_OK;}
    if(m==1&&!second&&fabs(phi)>=TE_HALF_PI)return TINYELLIPTIC_SINGULAR;
    if(fabs(phi)<=TE_HALF_PI)return principal(phi,m,second,value);
    n=floor((phi+TE_HALF_PI)/TE_PI);remainder=(phi-n*TE_PI);
    status=principal(remainder,m,second,&tail);if(status!=TINYELLIPTIC_OK)return status;
    if(n==0){*value=tail;return TINYELLIPTIC_OK;}
    status=principal(TE_HALF_PI,m,second,&complete);if(status!=TINYELLIPTIC_OK)return status;
    *value=2*n*complete+tail;return isfinite(*value)?TINYELLIPTIC_OK:TINYELLIPTIC_ARITHMETIC_ERROR;
}
static tinyelliptic_status evaluate(double phi,double m,int second,double *value){
    double out;tinyelliptic_status status;
    if(!value||!isfinite(phi)||!isfinite(m)||m<0||m>1)return TINYELLIPTIC_INVALID_ARGUMENT;
    if(fabs(phi)>1e6)return TINYELLIPTIC_RANGE;
    status=extended(phi,m,second,&out);if(status==TINYELLIPTIC_OK)*value=out;return status;
}
tinyelliptic_status tinyelliptic_first(double phi,double m,double *value){return evaluate(phi,m,0,value);}
tinyelliptic_status tinyelliptic_second(double phi,double m,double *value){return evaluate(phi,m,1,value);}
tinyelliptic_status tinyelliptic_complete_first(double m,double *value){return evaluate(TE_HALF_PI,m,0,value);}
tinyelliptic_status tinyelliptic_complete_second(double m,double *value){return evaluate(TE_HALF_PI,m,1,value);}
tinyelliptic_status tinyelliptic_inverse_second(double value,double m,double *amplitude){
    double complete,n,target,lo=-TE_HALF_PI,hi=TE_HALF_PI,x,result;
    if(!amplitude||!isfinite(value)||!isfinite(m)||m<0||m>1)return TINYELLIPTIC_INVALID_ARGUMENT;
    if(fabs(value)>1e6)return TINYELLIPTIC_RANGE;
    if(m==0){*amplitude=value;return TINYELLIPTIC_OK;}
    if(principal(TE_HALF_PI,m,1,&complete)!=TINYELLIPTIC_OK)return TINYELLIPTIC_NO_CONVERGENCE;
    if(fabs(value)<=complete){n=0;target=value;}
    else{n=floor((value+complete)/(2*complete));target=value-2*n*complete;}
    if(m==1)x=asin(fmax(-1,fmin(1,target)));
    else if(target<=-complete)x=-TE_HALF_PI;
    else if(target>=complete)x=TE_HALF_PI;
    else{
        int converged=0;x=target/complete*TE_HALF_PI;
        for(size_t iteration=0;iteration<80;++iteration){
            double e,derivative,candidate;
            if(principal(x,m,1,&e)!=TINYELLIPTIC_OK)return TINYELLIPTIC_NO_CONVERGENCE;
            e-=target;
            if(e==0){converged=1;break;}
            if(e>0)hi=x;else lo=x;
            derivative=sqrt(fmax(0,1-m*sin(x)*sin(x)));
            candidate=derivative>0?x-e/derivative:(lo+hi)/2;
            if(!isfinite(candidate)||candidate<=lo||candidate>=hi)candidate=(lo+hi)/2;
            if(fabs(candidate-x)<=8*DBL_EPSILON){x=candidate;converged=1;break;}
            x=candidate;
        }
        if(!converged)return TINYELLIPTIC_NO_CONVERGENCE;
    }
    result=n*TE_PI+x;if(!isfinite(result))return TINYELLIPTIC_ARITHMETIC_ERROR;
    if(fabs(result)>1e6)return TINYELLIPTIC_RANGE;
    *amplitude=result;return TINYELLIPTIC_OK;
}
tinyelliptic_status tinyelliptic_arc(double a,double b,double start,double end,double *length){
    double major,ratio,m,offset,n0,n1,r0,r1,e0,e1,complete,result;
    if(!length||!isfinite(a)||a<=0||!isfinite(b)||b<=0||!isfinite(start)||!isfinite(end))return TINYELLIPTIC_INVALID_ARGUMENT;
    if(fabs(start)>1e6||fabs(end)>1e6)return TINYELLIPTIC_RANGE;
    if(a==b){result=a*(end-start);if(!isfinite(result))return TINYELLIPTIC_ARITHMETIC_ERROR;*length=result;return TINYELLIPTIC_OK;}
    if(a>=b){major=a;ratio=b/a;offset=TE_HALF_PI;}else{major=b;ratio=a/b;offset=0;}
    m=1-ratio*ratio;n0=floor((start+offset+TE_HALF_PI)/TE_PI);n1=floor((end+offset+TE_HALF_PI)/TE_PI);
    r0=(start+offset)-n0*TE_PI;r1=(end+offset)-n1*TE_PI;
    if(principal(r0,m,1,&e0)!=TINYELLIPTIC_OK||principal(r1,m,1,&e1)!=TINYELLIPTIC_OK||principal(TE_HALF_PI,m,1,&complete)!=TINYELLIPTIC_OK)return TINYELLIPTIC_NO_CONVERGENCE;
    result=major*(2*(n1-n0)*complete+e1-e0);
    if(!isfinite(result))return TINYELLIPTIC_ARITHMETIC_ERROR;*length=result;return TINYELLIPTIC_OK;
}
