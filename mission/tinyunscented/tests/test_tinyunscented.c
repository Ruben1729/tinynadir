/* SPDX-License-Identifier: MIT */
#include "tinyunscented/tinyunscented.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "data/reference.h"
#define CHECK(x) do { if (!(x)) {fprintf(stderr,"line %d: %s\n",__LINE__,#x);return 1;} } while(0)
static int nonlinear(const double *x,size_t n,double *y,size_t m,void *context)
{size_t i;(void)context;if(n!=m)return 0;for(i=0;i<n;++i)y[i]=x[i]*x[i]+.1*x[(i+1)%n];return 1;}
int main(void)
{
    size_t i,j;double maximum=0;
    for(i=0;i<sizeof reference/sizeof reference[0];++i){size_t n=reference[i].n;tinyunscented_sigma sigma;double mean[6],P[36];CHECK(tinyunscented_generate(n,reference[i].mean,reference[i].P,reference[i].alpha,2,0,&sigma)==TINYUNSCENTED_OK);for(j=0;j<(2*n+1)*n;++j)CHECK(fabs(sigma.points[j]-reference[i].points[j])<2e-14);for(j=0;j<2*n+1;++j){CHECK(fabs(sigma.mean_weights[j]-reference[i].wm[j])<2e-14);CHECK(fabs(sigma.covariance_weights[j]-reference[i].wc[j])<2e-14);}CHECK(tinyunscented_statistics(sigma.count,n,sigma.points,sigma.mean_weights,sigma.covariance_weights,mean,P)==TINYUNSCENTED_OK);for(j=0;j<n;++j)CHECK(fabs(mean[j]-reference[i].mean[j])<2e-13);for(j=0;j<n*n;++j)CHECK(fabs(P[j]-reference[i].P[j])<3e-13);CHECK(tinyunscented_transform(&sigma,n,nonlinear,NULL,mean,P)==TINYUNSCENTED_OK);for(j=0;j<n;++j){double error=fabs(mean[j]-reference[i].out_mean[j]);if(error>maximum)maximum=error;CHECK(error<3e-12);}for(j=0;j<n*n;++j){double error=fabs(P[j]-reference[i].out_cov[j]);if(error>maximum)maximum=error;CHECK(error<3e-11);}}
    {double mean[2]={0,0},P[4]={1e20,0,0,1e-20},m[2],C[4];tinyunscented_sigma sigma,saved;CHECK(tinyunscented_generate(2,mean,P,.7,2,0,&sigma)==TINYUNSCENTED_OK);CHECK(tinyunscented_statistics(5,2,sigma.points,sigma.mean_weights,sigma.covariance_weights,m,C)==TINYUNSCENTED_OK);CHECK(fabs(C[0]/P[0]-1)<1e-14);CHECK(fabs(C[3]/P[3]-1)<1e-14);saved=sigma;P[3]=-1;CHECK(tinyunscented_generate(2,mean,P,.7,2,0,&sigma)!=TINYUNSCENTED_OK);CHECK(!memcmp(&sigma,&saved,sizeof sigma));P[0]=1;P[3]=1;P[1]=P[2]=2;CHECK(tinyunscented_generate(2,mean,P,.7,2,0,&sigma)==TINYUNSCENTED_NOT_PSD);}
    printf("Unscented: 9 FilterPy sigma/transform cases, dimensions2/3/6; max moment error %.9g\n",maximum);return 0;
}
