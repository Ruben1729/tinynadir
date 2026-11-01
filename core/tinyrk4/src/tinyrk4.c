/* SPDX-License-Identifier: MIT */
#include "tinyrk4/tinyrk4.h"
#include <math.h>
#include <string.h>
static int finite_values(const double *v,size_t n){if(!v)return 0;for(size_t i=0;i<n;++i)if(!isfinite(v[i]))return 0;return 1;}
static int evaluate(tinyrk4_rhs f,void *ctx,size_t n,double t,const double *y,double *k){
    if(!isfinite(t)||!finite_values(y,n))return 0;
    for(size_t i=0;i<n;++i)k[i]=NAN;
    return f(t,y,n,ctx,k)==0&&finite_values(k,n);
}
tinyrk4_status tinyrk4_step(tinyrk4_rhs f,void *ctx,size_t n,double t,double h,const double *y,double *out,double *w,size_t count){
    double *k1=w,*k2,*k3,*k4,*temp,*result;
    if(!f||!out||!w||n==0||n>TINYRK4_MAX_STATE||count<6*n||!isfinite(t)||!isfinite(h)||!finite_values(y,n))return TINYRK4_INVALID;
    if(!isfinite(t+h))return TINYRK4_ARITHMETIC;
    if(h==0){memmove(out,y,n*sizeof *out);return TINYRK4_OK;}
    if(t+h==t)return TINYRK4_ARITHMETIC;
    k2=w+n;k3=w+2*n;k4=w+3*n;temp=w+4*n;result=w+5*n;
    if(!evaluate(f,ctx,n,t,y,k1))return TINYRK4_DERIVATIVE_ERROR;
    for(size_t i=0;i<n;++i)temp[i]=y[i]+.5*h*k1[i];
    if(!evaluate(f,ctx,n,t+.5*h,temp,k2))return TINYRK4_DERIVATIVE_ERROR;
    for(size_t i=0;i<n;++i)temp[i]=y[i]+.5*h*k2[i];
    if(!evaluate(f,ctx,n,t+.5*h,temp,k3))return TINYRK4_DERIVATIVE_ERROR;
    for(size_t i=0;i<n;++i)temp[i]=y[i]+h*k3[i];
    if(!evaluate(f,ctx,n,t+h,temp,k4))return TINYRK4_DERIVATIVE_ERROR;
    for(size_t i=0;i<n;++i)result[i]=y[i]+h*(k1[i]/6+k2[i]/3+k3[i]/3+k4[i]/6);
    if(!finite_values(result,n))return TINYRK4_ARITHMETIC;
    memcpy(out,result,n*sizeof *out);return TINYRK4_OK;
}
