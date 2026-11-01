/* SPDX-License-Identifier: MIT */
#include "tinygravity/tinygravity.h"
#include <math.h>
#include <string.h>
/* Cartesian solid harmonics with forward-mode spatial derivatives avoid
   longitude/colatitude singularities at the poles. No spherical divisions. */
typedef struct {double v,d[3];} dual;
static dual constant(double x){dual r={x,{0,0,0}};return r;}
static dual add(dual a,dual b){a.v+=b.v;for(size_t i=0;i<3;++i)a.d[i]+=b.d[i];return a;}
static dual scale(dual a,double f){a.v*=f;for(size_t i=0;i<3;++i)a.d[i]*=f;return a;}
static dual multiply(dual a,dual b){dual r;r.v=a.v*b.v;for(size_t i=0;i<3;++i)r.d[i]=a.d[i]*b.v+a.v*b.d[i];return r;}
static dual reciprocal(dual a){dual r;r.v=1/a.v;for(size_t i=0;i<3;++i)r.d[i]=-a.d[i]*r.v*r.v;return r;}
static dual root(dual a){dual r;r.v=sqrt(a.v);for(size_t i=0;i<3;++i)r.d[i]=a.d[i]/(2*r.v);return r;}
static size_t index_nm(unsigned n,unsigned m){return (size_t)n*(n+1)/2+m;}
static void store(double *w,size_t k,dual a){w[4*k]=a.v;for(size_t j=0;j<3;++j)w[4*k+1+j]=a.d[j];}
static dual load(const double *w,size_t k){dual a;a.v=w[4*k];for(size_t j=0;j<3;++j)a.d[j]=w[4*k+1+j];return a;}
size_t tinygravity_coefficient_count(unsigned n){return n>TINYGRAVITY_MAX_DEGREE?0:((size_t)n+1)*(n+2)/2;}
size_t tinygravity_workspace_count(unsigned n){return 8*tinygravity_coefficient_count(n);}
/* Positive integer factorial ratio (n-m)!/(n+m)!. Degree<=32 bounds
 * the product to at most 64 divisions and keeps it above binary64 underflow.
 * Avoid lgamma: some platform implementations update shared signgam state. */
static double normalization(unsigned n,unsigned m){
    double ratio=1;
    for(unsigned k=n-m+1;k<=n+m;++k)ratio/=(double)k;
    return sqrt((m?2.0:1.0)*(2*n+1)*ratio);
}
tinygravity_status tinygravity_unnormalize(unsigned degree,const double *nc,const double *ns,size_t count,double *c,double *s){
    size_t required=tinygravity_coefficient_count(degree);
    if(!required||count<required||!nc||!ns||!c||!s||c==s)return TINYGRAVITY_INVALID;
    for(size_t i=0;i<required;++i)if(!isfinite(nc[i])||!isfinite(ns[i]))return TINYGRAVITY_INVALID;
    /* Validate all products before the first write, including in-place use. */
    for(unsigned n=0;n<=degree;++n)for(unsigned m=0;m<=n;++m){size_t k=index_nm(n,m);double norm=normalization(n,m);if(!isfinite(nc[k]*norm)||!isfinite(ns[k]*norm))return TINYGRAVITY_ARITHMETIC;}
    for(unsigned n=0;n<=degree;++n)for(unsigned m=0;m<=n;++m){size_t k=index_nm(n,m);double norm=normalization(n,m),cv=nc[k]*norm,sv=ns[k]*norm;c[k]=cv;s[k]=sv;}
    return TINYGRAVITY_OK;
}
tinygravity_status tinygravity_evaluate(double mu,double radius,unsigned degree,unsigned order,const double *c,const double *s,size_t count,const double r[3],double out[3],double *potential,double *work,size_t work_count){
    size_t required=tinygravity_coefficient_count(degree);double *vw=work,*ww,accel[3],p;dual x[3],r2,inverse,u,v,z,total;
    if(!required||order>degree||!r||!out||!c||!s||count<required||!work||work_count<8*required||!isfinite(mu)||mu<=0||!isfinite(radius)||radius<=0)return TINYGRAVITY_INVALID;
    for(size_t i=0;i<required;++i)if(!isfinite(c[i])||!isfinite(s[i]))return TINYGRAVITY_INVALID;
    r2=constant(0);for(size_t i=0;i<3;++i){if(!isfinite(r[i]))return TINYGRAVITY_INVALID;x[i]=constant(r[i]/radius);x[i].d[i]=1;r2=add(r2,multiply(x[i],x[i]));}
    if(r2.v==0)return TINYGRAVITY_SINGULAR;
    if(!isfinite(r2.v))return TINYGRAVITY_ARITHMETIC;
    inverse=reciprocal(r2);u=multiply(x[0],inverse);v=multiply(x[1],inverse);z=multiply(x[2],inverse);ww=work+4*required;
    store(vw,0,reciprocal(root(r2)));store(ww,0,constant(0));
    for(unsigned n=1;n<=degree;++n)for(unsigned m=0;m<=n&&m<=order;++m){dual a,b;
        if(m==n){dual prev_v=load(vw,index_nm(n-1,m-1)),prev_w=load(ww,index_nm(n-1,m-1));a=scale(add(multiply(u,prev_v),scale(multiply(v,prev_w),-1)),2*n-1);b=scale(add(multiply(u,prev_w),multiply(v,prev_v)),2*n-1);}
        else {double f=(2.0*n-1)/(n-m);a=scale(multiply(z,load(vw,index_nm(n-1,m))),f);b=scale(multiply(z,load(ww,index_nm(n-1,m))),f);if(n>=m+2){f=-(n+m-1.0)/(n-m);a=add(a,scale(multiply(inverse,load(vw,index_nm(n-2,m))),f));b=add(b,scale(multiply(inverse,load(ww,index_nm(n-2,m))),f));}}
        store(vw,index_nm(n,m),a);store(ww,index_nm(n,m),b);
    }
    total=constant(0);for(unsigned n=0;n<=degree;++n)for(unsigned m=0;m<=n&&m<=order;++m){size_t k=index_nm(n,m);total=add(total,add(scale(load(vw,k),c[k]),scale(load(ww,k),s[k])));}
    p=(mu/radius)*total.v;for(size_t i=0;i<3;++i){accel[i]=(mu/radius/radius)*total.d[i];if(!isfinite(accel[i]))return TINYGRAVITY_ARITHMETIC;}
    if(!isfinite(p))return TINYGRAVITY_ARITHMETIC;
    memcpy(out,accel,sizeof accel);if(potential)*potential=p;return TINYGRAVITY_OK;
}
