/* SPDX-License-Identifier: MIT */
#include "tinyhills/tinyhills.h"
#include <math.h>
#include <float.h>
static int matrix(double n,double dt,double P[36])
{
    double x,s,c,A,B,C,D;unsigned i;
    if(!isfinite(n) || n<0 || !isfinite(dt) || !isfinite(n*dt) || fabs(n*dt)>1e8) return 0;
    for(i=0;i<36;++i)P[i]=0;
    if(n==0) {
        for(i=0;i<6;++i)P[6*i+i]=1;
        for(i=0;i<3;++i)P[6*i+i+3]=dt;
        return 1;
    }
    x=n*dt;s=sin(x);c=cos(x);A=2*sin(x/2)*sin(x/2);
    B=fabs(x)<1e-5 ? dt*(1-x*x/6+x*x*x*x/120) : s/n;
    C=fabs(x)<1e-5 ? dt*(x/2-x*x*x/24+x*x*x*x*x/720) : A/n;
    if(fabs(x)<.1) {
        double x2=x*x;
        D=dt*x2*(-1.0/6+x2*(1.0/120+x2*(-1.0/5040+x2/362880)));
    } else D=(s-x)/n;
    P[0]=1+3*A;P[3]=B;P[4]=2*C;
    P[6]=6*n*D;P[7]=1;P[9]=-2*C;P[10]=dt+4*D;
    P[14]=c;P[17]=B;
    P[18]=3*n*s;P[21]=c;P[22]=2*s;
    P[24]=-6*n*A;P[27]=-2*s;P[28]=1-4*A;
    P[32]=-n*s;P[35]=c;
    for(i=0;i<36;++i)if(!isfinite(P[i]))return 0;
    return 1;
}
tinyhills_status tinyhills_transition(double n,double dt,double out[36])
{
    double P[36];unsigned i;
    if(!out || !matrix(n,dt,P))return TINYHILLS_INVALID;
    for(i=0;i<36;++i)out[i]=P[i];
    return TINYHILLS_OK;
}
tinyhills_status tinyhills_propagate(double n,double dt,const tinyhills_state *initial,tinyhills_state *out)
{
    double P[36],input[6],result[6];tinyhills_state state;unsigned i,j;
    if(!initial || !out || !matrix(n,dt,P))return TINYHILLS_INVALID;
    for(i=0;i<3;++i){input[i]=initial->position_m[i];input[3+i]=initial->velocity_m_s[i];}
    for(i=0;i<6;++i)if(!isfinite(input[i]))return TINYHILLS_INVALID;
    for(i=0;i<6;++i){result[i]=0;for(j=0;j<6;++j)result[i]+=P[6*i+j]*input[j];if(!isfinite(result[i]))return TINYHILLS_INVALID;}
    for(i=0;i<3;++i){state.position_m[i]=result[i];state.velocity_m_s[i]=result[3+i];}
    *out=state;return TINYHILLS_OK;
}
tinyhills_status tinyhills_rendezvous(double n,double dt,const double r[3],double out[3])
{
    double P[36],a,b,c,d,scale,det,y0,y1,y2,v[3];unsigned j;
    if(!r || !out || dt==0 || !matrix(n,dt,P))return TINYHILLS_INVALID;
    for(j=0;j<3;++j)if(!isfinite(r[j]))return TINYHILLS_INVALID;
    a=P[3];b=P[4];c=P[9];d=P[10];scale=fmax(fmax(fabs(a),fabs(b)),fmax(fabs(c),fabs(d)));
    if(scale==0 || fabs(P[17])<=64*DBL_EPSILON*fabs(dt))return TINYHILLS_SINGULAR;
    a/=scale;b/=scale;c/=scale;d/=scale;det=a*d-b*c;
    if(fabs(det)<1e-12)return TINYHILLS_SINGULAR;
    y0=-(P[0]*r[0]+P[1]*r[1]+P[2]*r[2])/scale;
    y1=-(P[6]*r[0]+P[7]*r[1]+P[8]*r[2])/scale;
    y2=-(P[12]*r[0]+P[13]*r[1]+P[14]*r[2]);
    v[0]=(d*y0-b*y1)/det;v[1]=(-c*y0+a*y1)/det;v[2]=y2/P[17];
    for(j=0;j<3;++j)if(!isfinite(v[j]))return TINYHILLS_INVALID;
    for(j=0;j<3;++j)out[j]=v[j];
    return TINYHILLS_OK;
}
