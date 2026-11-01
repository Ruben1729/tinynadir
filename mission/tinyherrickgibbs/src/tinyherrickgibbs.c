/* SPDX-License-Identifier: MIT */
#include "tinyherrickgibbs/tinyherrickgibbs.h"
#include <math.h>
static double norm(const double x[3]){return hypot(hypot(x[0],x[1]),x[2]);}
static double dot(const double a[3],const double b[3]){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
static void cross(const double a[3],const double b[3],double c[3]){c[0]=a[1]*b[2]-a[2]*b[1];c[1]=a[2]*b[0]-a[0]*b[2];c[2]=a[0]*b[1]-a[1]*b[0];}
tinyherrickgibbs_status tinyherrickgibbs_estimate(double mu,const double r1[3],const double r2[3],const double r3[3],double a,double b,double angle_tol,double plane_tol,tinyherrickgibbs_result *out)
{
    const double *src[3]={r1,r2,r3};double r[3][3],length[3],scale,c12[3],c23[3],cn,plane,accel,k[3];
    tinyherrickgibbs_result result;unsigned i,j;
    if(!r1 || !r2 || !r3 || !out || !isfinite(mu) || mu<=0 || !isfinite(a) || a<=0 || !isfinite(b) || b<=0 ||
       !isfinite(a+b) || !isfinite(a*b) || !isfinite(a*(a+b)) || !isfinite(b*(a+b)) ||
       !isfinite(angle_tol) || angle_tol<=0 || angle_tol>1.5707963267948966 ||
       !isfinite(plane_tol) || plane_tol<=0 || plane_tol>1.5707963267948966)return TINYHERRICKGIBBS_INVALID;
    for(i=0;i<3;++i){for(j=0;j<3;++j)if(!isfinite(src[i][j]))return TINYHERRICKGIBBS_INVALID;length[i]=norm(src[i]);if(!(length[i]>0) || !isfinite(length[i]))return TINYHERRICKGIBBS_INVALID;}
    scale=fmax(length[0],fmax(length[1],length[2]));
    for(i=0;i<3;++i){length[i]/=scale;for(j=0;j<3;++j)r[i][j]=src[i][j]/scale;}
    cross(r[0],r[1],c12);cross(r[1],r[2],c23);cn=norm(c23);
    if(cn<1e-14 || norm(c12)<1e-14)return TINYHERRICKGIBBS_GEOMETRY;
    plane=dot(r[0],c23)/(length[0]*cn);plane=fmax(-1.0,fmin(1.0,plane));result.coplanarity_rad=asin(plane);
    result.angle12_rad=atan2(norm(c12),dot(r[0],r[1]));result.angle23_rad=atan2(cn,dot(r[1],r[2]));
    if(fabs(result.coplanarity_rad)>plane_tol || result.angle12_rad>angle_tol || result.angle23_rad>angle_tol)return TINYHERRICKGIBBS_GEOMETRY;
    accel=(mu/scale)/scale;
    k[0]=-b/(a*(a+b))*scale-accel*b/(12*length[0]*length[0]*length[0]);
    k[1]=(b-a)/(a*b)*scale+accel*(b-a)/(12*length[1]*length[1]*length[1]);
    k[2]=a/(b*(a+b))*scale+accel*a/(12*length[2]*length[2]*length[2]);
    for(j=0;j<3;++j){result.velocity_m_s[j]=k[0]*r[0][j]+k[1]*r[1][j]+k[2]*r[2][j];if(!isfinite(result.velocity_m_s[j]))return TINYHERRICKGIBBS_INVALID;}
    *out=result;return TINYHERRICKGIBBS_OK;
}
