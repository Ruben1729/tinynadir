/* SPDX-License-Identifier: MIT
   Independent spherical-harmonic synthesis; official IAGA coefficients. */
#include "tinyigrf/tinyigrf.h"
#include <math.h>
#include <stddef.h>
#include "coefficients.h"
#define PI 3.141592653589793238462643383279502884

tinyigrf_status tinyigrf_geocentric(double year,double radius,double theta,double phi,unsigned degree,double out[3]) {
    double p[14][14]={{0}},dp[14][14]={{0}},cl[14],sl[14];
    double s,c,ar,power,br=0,bt=0,bp=0;unsigned n,m;
    if(!out || !isfinite(year)||!isfinite(radius)||!isfinite(theta)||!isfinite(phi))return TINYIGRF_INVALID;
    if(year<2025 || year>2030 || radius<6000000 || theta<0 || theta>PI || degree<1 || degree>13)return TINYIGRF_RANGE;
    phi=remainder(phi,2*PI);s=sin(theta);c=cos(theta);ar=6371200.0/radius;power=ar*ar*ar;
    cl[0]=1;sl[0]=0;
    for(m=1;m<=degree;m++){cl[m]=cos(m*phi);sl[m]=sin(m*phi);}
    p[0][0]=1;
    for(n=1;n<=degree;n++) {
        for(m=0;m<=n;m++) {
            double norm=1,over,comb,g,h;unsigned j,index=n*(n+1)/2+m;
            if(m==n) {
                p[n][m]=(2.0*n-1)*s*p[n-1][m-1];
                dp[n][m]=(2.0*n-1)*(c*p[n-1][m-1]+s*dp[n-1][m-1]);
            } else if(n==m+1) {
                p[n][m]=(2.0*n-1)*c*p[n-1][m];
                dp[n][m]=(2.0*n-1)*(-s*p[n-1][m]+c*dp[n-1][m]);
            } else {
                p[n][m]=((2.0*n-1)*c*p[n-1][m]-(n+m-1.0)*p[n-2][m])/(n-m);
                dp[n][m]=((2.0*n-1)*(-s*p[n-1][m]+c*dp[n-1][m])-(n+m-1.0)*dp[n-2][m])/(n-m);
            }
            if(m) {for(j=n-m+1;j<=n+m;j++)norm/=j;norm=sqrt(2*norm);}
            g=igrf14_2025[index].g+(year-2025)*igrf14_2025[index].dg;
            h=igrf14_2025[index].h+(year-2025)*igrf14_2025[index].dh;
            comb=g*cl[m]+h*sl[m];
            br+=(n+1)*power*comb*p[n][m]*norm;
            bt-=power*comb*dp[n][m]*norm;
            if(m) {
                if(fabs(s)>1e-12)over=p[n][m]/s;
                else if(m==1)over=.5*n*(n+1)*((c>=0 || n%2)?1:-1);
                else over=0;
                bp+=power*m*(g*sl[m]-h*cl[m])*over*norm;
            }
        }
        power*=ar;
    }
    if(!isfinite(br)||!isfinite(bt)||!isfinite(bp))return TINYIGRF_RANGE;
    out[0]=br*1e-9;out[1]=bt*1e-9;out[2]=bp*1e-9;return TINYIGRF_OK;
}
tinyigrf_status tinyigrf_ecef(double year,const double r[3],unsigned degree,double out[3]) {
    double v[3],b[3],radius,theta,phi,s,c,sp,cp; tinyigrf_status status;
    if(!r||!out||!isfinite(r[0])||!isfinite(r[1])||!isfinite(r[2]))return TINYIGRF_INVALID;
    radius=hypot(hypot(r[0],r[1]),r[2]);
    if(!isfinite(radius)||radius<6000000)return TINYIGRF_RANGE;
    theta=atan2(hypot(r[0],r[1]),r[2]);phi=atan2(r[1],r[0]);
    status=tinyigrf_geocentric(year,radius,theta,phi,degree,b);if(status!=TINYIGRF_OK)return status;
    s=sin(theta);c=cos(theta);sp=sin(phi);cp=cos(phi);
    v[0]=b[0]*s*cp+b[1]*c*cp-b[2]*sp;
    v[1]=b[0]*s*sp+b[1]*c*sp+b[2]*cp;
    v[2]=b[0]*c-b[1]*s;
    out[0]=v[0];out[1]=v[1];out[2]=v[2];return TINYIGRF_OK;
}
tinyigrf_status tinyigrf_geodetic(double year,double lat,double lon,double height,unsigned degree,double out[3]) {
    double r[3],s,c,n;const double a=6378137.0,e2=0.0066943799901413165;
    if(!out||!isfinite(lat)||!isfinite(lon)||!isfinite(height))return TINYIGRF_INVALID;
    if(lat < -PI/2 || lat > PI/2 || height < -1000)return TINYIGRF_RANGE;
    lon=remainder(lon,2*PI);s=sin(lat);c=cos(lat);n=a/sqrt(1-e2*s*s);
    r[0]=(n+height)*c*cos(lon);r[1]=(n+height)*c*sin(lon);r[2]=(n*(1-e2)+height)*s;
    return tinyigrf_ecef(year,r,degree,out);
}
