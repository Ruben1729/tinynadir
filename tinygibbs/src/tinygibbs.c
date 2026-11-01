/* SPDX-License-Identifier: MIT */
#include "tinygibbs/tinygibbs.h"
#include <math.h>
static double norm(const double x[3]){return hypot(hypot(x[0],x[1]),x[2]);}
static double dot(const double a[3],const double b[3]){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
static void cross(const double a[3],const double b[3],double c[3]){c[0]=a[1]*b[2]-a[2]*b[1];c[1]=a[2]*b[0]-a[0]*b[2];c[2]=a[0]*b[1]-a[1]*b[0];}
tinygibbs_status tinygibbs_estimate(double mu,const double r1[3],const double r2[3],const double r3[3],double tol,tinygibbs_result *out)
{
    double r[3][3],length[3],scale,c12[3],c23[3],c31[3],N[3],D[3],S[3],B[3],nn,dn,cn,factor,plane;
    tinygibbs_result result;const double *src[3]={r1,r2,r3};unsigned i,j;
    if(!r1 || !r2 || !r3 || !out || !isfinite(mu) || mu<=0 || !isfinite(tol) || tol<0 || tol>1.5707963267948966)return TINYGIBBS_INVALID;
    for(i=0;i<3;++i){for(j=0;j<3;++j)if(!isfinite(src[i][j]))return TINYGIBBS_INVALID;length[i]=norm(src[i]);if(!(length[i]>0) || !isfinite(length[i]))return TINYGIBBS_INVALID;}
    scale=fmax(length[0],fmax(length[1],length[2]));
    for(i=0;i<3;++i){length[i]/=scale;for(j=0;j<3;++j)r[i][j]=src[i][j]/scale;}
    cross(r[0],r[1],c12);cross(r[1],r[2],c23);cross(r[2],r[0],c31);cn=norm(c23);
    if(cn<1e-14 || norm(c12)<1e-14)return TINYGIBBS_GEOMETRY;
    plane=dot(r[0],c23)/(length[0]*cn);plane=fmax(-1.0,fmin(1.0,plane));
    result.coplanarity_rad=asin(plane);
    if(fabs(result.coplanarity_rad)>tol)return TINYGIBBS_GEOMETRY;
    result.angle12_rad=atan2(norm(c12),dot(r[0],r[1]));result.angle23_rad=atan2(cn,dot(r[1],r[2]));
    for(j=0;j<3;++j){N[j]=length[0]*c23[j]+length[1]*c31[j]+length[2]*c12[j];D[j]=c12[j]+c23[j]+c31[j];S[j]=(length[1]-length[2])*r[0][j]+(length[2]-length[0])*r[1][j]+(length[0]-length[1])*r[2][j];}
    nn=norm(N);dn=norm(D);
    if(nn<1e-14 || dn<1e-14 || dot(N,D)<=0)return TINYGIBBS_GEOMETRY;
    cross(D,r[1],B);factor=sqrt(mu/scale)/sqrt(nn*dn);
    for(j=0;j<3;++j){result.velocity_m_s[j]=factor*(B[j]/length[1]+S[j]);if(!isfinite(result.velocity_m_s[j]))return TINYGIBBS_INVALID;}
    *out=result;return TINYGIBBS_OK;
}
