/* SPDX-License-Identifier: MIT */
#include "tinyshadow/tinyshadow.h"
#include <math.h>
#define TS_PI 3.14159265358979323846264338327950288
static int finite3(const double *v){return v&&isfinite(v[0])&&isfinite(v[1])&&isfinite(v[2]);}
static double norm3(const double v[3]){return hypot(hypot(v[0],v[1]),v[2]);}
static double clamp(double x){return fmax(-1,fmin(1,x));}
tinyshadow_status tinyshadow_sun_fraction(const double r[3],const double sun[3],double rb,double rs,double *out){
    double s[3],earth[3],cross[3],rnorm,d,a,b,delta,dot,overlap,result;
    if(!out||!finite3(r)||!finite3(sun)||!isfinite(rb)||rb<=0||!isfinite(rs)||rs<=0)return TINYSHADOW_INVALID;
    for(size_t i=0;i<3;++i)s[i]=sun[i]-r[i];
    rnorm=norm3(r);d=norm3(s);
    if(!isfinite(rnorm)||!isfinite(d))return TINYSHADOW_ARITHMETIC;
    if(rnorm<=rb||d<=rs)return TINYSHADOW_GEOMETRY;
    for(size_t i=0;i<3;++i){s[i]/=d;earth[i]=-r[i]/rnorm;}
    cross[0]=s[1]*earth[2]-s[2]*earth[1];cross[1]=s[2]*earth[0]-s[0]*earth[2];cross[2]=s[0]*earth[1]-s[1]*earth[0];
    dot=s[0]*earth[0]+s[1]*earth[1]+s[2]*earth[2];delta=atan2(norm3(cross),dot);a=asin(rs/d);b=asin(rb/rnorm);
    if(delta>=a+b)result=1;
    else if(b>=delta+a)result=0;
    else if(a>=delta+b)result=1-(b/a)*(b/a);
    else {double x=(delta*delta+a*a-b*b)/(2*delta*a),y=(delta*delta+b*b-a*a)/(2*delta*b);double product=(-delta+a+b)*(delta+a-b)*(delta-a+b)*(delta+a+b);overlap=a*a*acos(clamp(x))+b*b*acos(clamp(y))-.5*sqrt(fmax(0,product));result=1-overlap/(TS_PI*a*a);}
    if(!isfinite(result))return TINYSHADOW_ARITHMETIC;
    *out=fmax(0,fmin(1,result));return TINYSHADOW_OK;
}
tinyshadow_status tinyshadow_cylindrical(const double r[3],const double sun[3],double radius,int *out){
    double n,dot,perp[3];if(!out||!finite3(r)||!finite3(sun)||!isfinite(radius)||radius<=0)return TINYSHADOW_INVALID;
    n=norm3(sun);if(n==0||!isfinite(n))return TINYSHADOW_GEOMETRY;
    dot=r[0]*(sun[0]/n)+r[1]*(sun[1]/n)+r[2]*(sun[2]/n);for(size_t i=0;i<3;++i)perp[i]=r[i]-dot*(sun[i]/n);
    if(!isfinite(dot)||!isfinite(norm3(perp)))return TINYSHADOW_ARITHMETIC;
    *out=dot<0&&norm3(perp)<=radius;return TINYSHADOW_OK;
}
tinyshadow_status tinyshadow_segment_ellipsoid(const double from[3],const double to[3],const double axes[3],int *out){
    double a[3],d[3],projection=0,denom=0,fraction,distance2=0;
    if(!out||!finite3(from)||!finite3(to)||!finite3(axes))return TINYSHADOW_INVALID;
    for(size_t i=0;i<3;++i){if(axes[i]<=0)return TINYSHADOW_INVALID;a[i]=from[i]/axes[i];d[i]=to[i]/axes[i]-a[i];projection+=a[i]*d[i];denom+=d[i]*d[i];}
    if(!isfinite(denom)||!isfinite(projection))return TINYSHADOW_ARITHMETIC;
    fraction=denom?fmax(0,fmin(1,-projection/denom)):0;for(size_t i=0;i<3;++i){double v=a[i]+fraction*d[i];distance2+=v*v;}
    if(!isfinite(distance2))return TINYSHADOW_ARITHMETIC;
    *out=distance2<=1;return TINYSHADOW_OK;
}
static int callback_value(tinyshadow_fraction_callback f,void *ctx,double t,double *value){return f(t,ctx,value)==0&&isfinite(*value)&&*value>=0&&*value<=1;}
tinyshadow_status tinyshadow_event(tinyshadow_fraction_callback f,void *ctx,double lo,double hi,double threshold,double tolerance,unsigned max_iterations,double *out){
    double a,b;
    if(!f||!out||!isfinite(lo)||!isfinite(hi)||hi<lo||!isfinite(hi-lo)||!isfinite(threshold)||threshold<0||threshold>1||!isfinite(tolerance)||tolerance<=0||!max_iterations||max_iterations>256)return TINYSHADOW_INVALID;
    if(!callback_value(f,ctx,lo,&a)||!callback_value(f,ctx,hi,&b))return TINYSHADOW_CALLBACK_ERROR;
    a-=threshold;b-=threshold;
    if(a==0){*out=lo;return TINYSHADOW_OK;}if(b==0){*out=hi;return TINYSHADOW_OK;}
    if((a>0)==(b>0))return TINYSHADOW_UNBRACKETED;
    for(unsigned i=0;i<max_iterations;++i){double mid=lo+.5*(hi-lo),value;if(hi-lo<=tolerance){*out=mid;return TINYSHADOW_OK;}if(mid==lo||mid==hi)return TINYSHADOW_NO_CONVERGENCE;if(!callback_value(f,ctx,mid,&value))return TINYSHADOW_CALLBACK_ERROR;value-=threshold;if(value==0){*out=mid;return TINYSHADOW_OK;}if((value>0)==(a>0)){lo=mid;a=value;}else hi=mid;}
    return TINYSHADOW_NO_CONVERGENCE;
}
