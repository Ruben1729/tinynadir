/* SPDX-License-Identifier: MIT */
#include "tinyinterpolation/tinyinterpolation.h"
#include <math.h>
static int valid(const double *x,const double *y,size_t n,size_t limit,double q,double *out) {
    if(!x||!y||!out||n<2||n>limit||!isfinite(q))return 0;
    for(size_t i=0;i<n;++i)if(!isfinite(x[i])||!isfinite(y[i])||(i&&x[i]<=x[i-1]))return 0;
    return 1;
}
tinyinterpolation_status tinyinterpolation_linear(const double *x,const double *y,size_t n,double q,double *out) {
    size_t lo=0,hi;double f,r,span;
    if(!valid(x,y,n,4096,q,out))return TINYINTERPOLATION_INVALID;
    if(q<x[0]||q>x[n-1])return TINYINTERPOLATION_RANGE;
    hi=n-1;
    while(hi-lo>1){size_t m=lo+(hi-lo)/2;if(q<x[m])hi=m;else lo=m;}
    span=x[hi]-x[lo];f=(q-x[lo])/span;r=(1-f)*y[lo]+f*y[hi];
    if(!isfinite(span)||!isfinite(r))return TINYINTERPOLATION_ARITHMETIC;
    *out=r;return TINYINTERPOLATION_OK;
}
tinyinterpolation_status tinyinterpolation_lagrange(const double *x,const double *y,size_t n,double q,double *out) {
    double scale,r=0,coordinate[8],z;
    if(!valid(x,y,n,8,q,out))return TINYINTERPOLATION_INVALID;
    if(q<x[0]||q>x[n-1])return TINYINTERPOLATION_RANGE;
    scale=x[n-1]-x[0];if(!isfinite(scale))return TINYINTERPOLATION_ARITHMETIC;
    for(size_t i=0;i<n;++i){if(q==x[i]){*out=y[i];return TINYINTERPOLATION_OK;}coordinate[i]=(x[i]-x[0])/scale;}
    z=(q-x[0])/scale;
    for(size_t i=0;i<n;++i){double term=y[i];for(size_t j=0;j<n;++j)if(i!=j)term*=(z-coordinate[j])/(coordinate[i]-coordinate[j]);r+=term;}
    if(!isfinite(r))return TINYINTERPOLATION_ARITHMETIC;
    *out=r;return TINYINTERPOLATION_OK;
}
tinyinterpolation_status tinyinterpolation_hermite(double x0,double x1,double y0,double y1,double m0,double m1,double q,double *out,double *deriv) {
    double h,t,t2,t3,r,d;
    if(!out||!isfinite(x0)||!isfinite(x1)||x1<=x0||!isfinite(y0)||!isfinite(y1)||!isfinite(m0)||!isfinite(m1)||!isfinite(q))return TINYINTERPOLATION_INVALID;
    if(q<x0||q>x1)return TINYINTERPOLATION_RANGE;
    h=x1-x0;t=(q-x0)/h;t2=t*t;t3=t2*t;
    r=(2*t3-3*t2+1)*y0+(t3-2*t2+t)*h*m0+(-2*t3+3*t2)*y1+(t3-t2)*h*m1;
    d=((6*t2-6*t)*y0+(-6*t2+6*t)*y1)/h+(3*t2-4*t+1)*m0+(3*t2-2*t)*m1;
    if(!isfinite(h)||!isfinite(r)||!isfinite(d))return TINYINTERPOLATION_ARITHMETIC;
    *out=r;if(deriv)*deriv=d;return TINYINTERPOLATION_OK;
}
