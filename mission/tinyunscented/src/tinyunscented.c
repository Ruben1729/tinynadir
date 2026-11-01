/* SPDX-License-Identifier: MIT */
#include "tinyunscented/tinyunscented.h"
#include <math.h>
#include <float.h>
static void add(double value,double *sum,double *correction)
{double y=value-*correction,next=*sum+y;*correction=(next-*sum)-y;*sum=next;}
tinyunscented_status tinyunscented_generate(size_t n,const double *mean,const double *P,double alpha,double beta,double kappa,tinyunscented_sigma *out)
{
    double L[36]={0},diagonal[6],c,tol;size_t i,j,k;tinyunscented_sigma result={0};
    if(!mean||!P||!out||n<1||n>TINYUNSCENTED_MAX_DIM||!isfinite(alpha)||alpha<=0||!isfinite(beta)||beta<0||!isfinite(kappa)||!((double)n+kappa>0))return TINYUNSCENTED_INVALID;
    c=alpha*alpha*((double)n+kappa);if(!isfinite(c)||!(c>0))return TINYUNSCENTED_INVALID;
    tol=64*DBL_EPSILON;
    for(i=0;i<n;++i){
        if(!isfinite(mean[i]) || !isfinite(P[n*i+i]) || P[n*i+i]<0)return TINYUNSCENTED_INVALID;
        diagonal[i]=sqrt(P[n*i+i]);
        for(j=0;j<n;++j)if(!isfinite(P[n*i+j]))return TINYUNSCENTED_INVALID;
    }
    for(i=0;i<n;++i)for(j=0;j<=i;++j){
        double value,reverse;
        if(diagonal[i]==0 || diagonal[j]==0){
            if(P[n*i+j]!=0 || P[n*j+i]!=0)return TINYUNSCENTED_NOT_PSD;
            value=0;
        } else {
            value=(P[n*i+j]/diagonal[i])/diagonal[j];
            reverse=(P[n*j+i]/diagonal[j])/diagonal[i];
            if(!isfinite(value) || !isfinite(reverse) || fabs(value-reverse)>tol*fmax(1,fmax(fabs(value),fabs(reverse))))return TINYUNSCENTED_INVALID;
            value=.5*value+.5*reverse;
            if(fabs(value)>1+tol)return TINYUNSCENTED_NOT_PSD;
        }
        for(k=0;k<j;++k)value-=L[n*i+k]*L[n*j+k];
        if(i==j){if(value< -tol)return TINYUNSCENTED_NOT_PSD;L[n*i+j]=value<=0 ? 0 : sqrt(value);}
        else if(L[n*j+j]==0){if(fabs(value)>tol)return TINYUNSCENTED_NOT_PSD;L[n*i+j]=0;}
        else L[n*i+j]=value/L[n*j+j];
    }
    result.dimension=n;result.count=2*n+1;
    for(i=0;i<result.count;++i){result.mean_weights[i]=0.5/c;result.covariance_weights[i]=0.5/c;}
    result.mean_weights[0]=1-(double)n/c;result.covariance_weights[0]=result.mean_weights[0]+1-alpha*alpha+beta;
    for(j=0;j<n;++j){result.points[j]=mean[j];for(i=0;i<n;++i){double offset=(sqrt(c)*diagonal[j])*L[n*j+i];result.points[(i+1)*n+j]=mean[j]+offset;result.points[(i+1+n)*n+j]=mean[j]-offset;if(offset!=0 && (result.points[(i+1)*n+j]==mean[j] || result.points[(i+1+n)*n+j]==mean[j]))return TINYUNSCENTED_INVALID;}}
    for(i=0;i<result.count;++i){if(!isfinite(result.mean_weights[i])||!isfinite(result.covariance_weights[i]))return TINYUNSCENTED_INVALID;for(j=0;j<n;++j)if(!isfinite(result.points[i*n+j]))return TINYUNSCENTED_INVALID;}
    *out=result;return TINYUNSCENTED_OK;
}
tinyunscented_status tinyunscented_statistics(size_t count,size_t n,const double *points,const double *wm,const double *wc,double *mean,double *cov)
{
    double m[6]={0},mc[6]={0},P[36]={0},pc[36]={0},wsum=0,wc_sum=0;size_t i,j,k;
    if(!points||!wm||!wc||!mean||!cov||count<1||count>TINYUNSCENTED_MAX_POINTS||n<1||n>TINYUNSCENTED_MAX_DIM)return TINYUNSCENTED_INVALID;
    for(i=0;i<count;++i){if(!isfinite(wm[i])||!isfinite(wc[i]))return TINYUNSCENTED_INVALID;add(wm[i],&wsum,&wc_sum);for(j=0;j<n;++j){double value=wm[i]*points[i*n+j];if(!isfinite(points[i*n+j])||!isfinite(value))return TINYUNSCENTED_INVALID;add(value,&m[j],&mc[j]);}}
    if(fabs(wsum-1)>1e-8)return TINYUNSCENTED_INVALID;
    for(i=0;i<count;++i)for(j=0;j<n;++j)for(k=0;k<=j;++k){double value=wc[i]*(points[i*n+j]-m[j])*(points[i*n+k]-m[k]);if(!isfinite(value))return TINYUNSCENTED_INVALID;add(value,&P[j*n+k],&pc[j*n+k]);}
    for(j=0;j<n;++j){if(!isfinite(m[j]))return TINYUNSCENTED_INVALID;for(k=0;k<=j;++k){if(!isfinite(P[j*n+k]))return TINYUNSCENTED_INVALID;P[k*n+j]=P[j*n+k];}}
    for(j=0;j<n;++j)mean[j]=m[j];
    for(j=0;j<n*n;++j)cov[j]=P[j];
    return TINYUNSCENTED_OK;
}
tinyunscented_status tinyunscented_transform(const tinyunscented_sigma *sigma,size_t m,tinyunscented_model model,void *context,double *mean,double *cov)
{
    double transformed[78],result_mean[6],result_cov[36];size_t i,j;tinyunscented_status status;
    if(!sigma||!model||!mean||!cov||sigma->dimension<1||sigma->dimension>6||sigma->count!=2*sigma->dimension+1||m<1||m>6)return TINYUNSCENTED_INVALID;
    for(i=0;i<sigma->count;++i){for(j=0;j<sigma->dimension;++j)if(!isfinite(sigma->points[i*sigma->dimension+j]))return TINYUNSCENTED_INVALID;for(j=0;j<m;++j)transformed[i*m+j]=NAN;if(!model(sigma->points+i*sigma->dimension,sigma->dimension,transformed+i*m,m,context))return TINYUNSCENTED_MODEL_ERROR;for(j=0;j<m;++j)if(!isfinite(transformed[i*m+j]))return TINYUNSCENTED_MODEL_ERROR;}
    status=tinyunscented_statistics(sigma->count,m,transformed,sigma->mean_weights,sigma->covariance_weights,result_mean,result_cov);
    if(status!=TINYUNSCENTED_OK)return status;
    for(i=0;i<m;++i)mean[i]=result_mean[i];
    for(i=0;i<m*m;++i)cov[i]=result_cov[i];
    return TINYUNSCENTED_OK;
}
