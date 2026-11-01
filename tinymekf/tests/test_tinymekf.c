#include "tinymekf/tinymekf.h"
#include "tinylinalg/tinylinalg.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL %d: %s\n",__LINE__,#c);return 1;}}while(0)
#include "data/reference.h"
static double distance(const double a[4],const double b[4]){double m=0,p=0;for(size_t i=0;i<4;++i){m+=(a[i]-b[i])*(a[i]-b[i]);p+=(a[i]+b[i])*(a[i]+b[i]);}return sqrt(fmin(m,p));}
static int covariance_same(const double a[36],const double b[36]){for(size_t i=0;i<36;++i)if(fabs(a[i]-b[i])>5e-17+3e-11*fabs(b[i]))return 0;return 1;}
int main(void){tinymekf_workspace work;tinymekf_state state;double nis,max_q=0,max_covariance=0;
 CHECK(tinymekf_initialize(&state,initial_q,initial_bias,initial_covariance,&work)==TINYMEKF_OK);
 for(size_t c=0;c<sizeof reference_cases/sizeof reference_cases[0];++c){
  CHECK(tinymekf_predict(&state,reference_cases[c].gyro,.05,gyro_psd,bias_psd,&work)==TINYMEKF_OK);
  CHECK(distance(state.q,reference_cases[c].predicted_q)<3e-11);CHECK(covariance_same(state.covariance,reference_cases[c].predicted_covariance));
  for(size_t k=0;k<2;++k){CHECK(tinymekf_update_vector(&state,reference_directions+3*k,reference_cases[c].body+3*k,1e-5,0,&nis,&work)==TINYMEKF_OK);double error=distance(state.q,reference_cases[c].update[k].q);CHECK(error<3e-11);max_q=fmax(max_q,error);for(size_t i=0;i<3;++i)CHECK(fabs(state.bias[i]-reference_cases[c].update[k].bias[i])<3e-12);CHECK(covariance_same(state.covariance,reference_cases[c].update[k].covariance));for(size_t i=0;i<36;++i)max_covariance=fmax(max_covariance,fabs(state.covariance[i]-reference_cases[c].update[k].covariance[i]));CHECK(fabs(nis-reference_cases[c].update[k].nis)<3e-10*fmax(1,nis));double lower[36],scratch[36];CHECK(tinylinalg_cholesky(6,state.covariance,lower,scratch,36)==TINYLINALG_OK);}
 }
 CHECK(distance(state.q,final_truth)<.01);
 {tinymekf_state before=state;double ref[]={1,0,0},opposite[]={-1,0,0};nis=123;CHECK(tinymekf_update_vector(&state,ref,opposite,1e-8,.001,&nis,&work)==TINYMEKF_REJECTED);CHECK(!memcmp(&before,&state,sizeof state)&&nis==123);CHECK(tinymekf_predict(&state,ref,-1,gyro_psd,bias_psd,&work)==TINYMEKF_INVALID_ARGUMENT);CHECK(!memcmp(&before,&state,sizeof state));}
 {double q[]={1,0,0,0},zero[]={0,0,0},p[36]={0},ref[]={0,0,1};for(size_t i=0;i<6;++i)p[i*6+i]=i<3?.01:.001;CHECK(tinymekf_initialize(&state,q,zero,p,&work)==TINYMEKF_OK);CHECK(tinymekf_update_vector(&state,ref,ref,.001,0,NULL,&work)==TINYMEKF_OK);CHECK(state.covariance[2*6+2]==.01);CHECK(state.covariance[0]<.01&&state.covariance[7]<.01);CHECK(distance(state.q,q)<1e-15);}
 {double q[]={1,0,0,0},zero[]={0,0,0},p[36]={0};for(size_t i=0;i<6;++i)p[i*6+i]=i<3?.01:.001;CHECK(tinymekf_initialize(&state,q,zero,p,&work)==TINYMEKF_OK);CHECK(tinymekf_predict(&state,zero,.1,zero,zero,&work)==TINYMEKF_OK);CHECK(fabs(state.covariance[0]-(.01+.001*.01))<1e-16);CHECK(fabs(state.covariance[3]+.0001)<1e-16);}
 {double q[]={1,0,0,0},zero[]={0,0,0},bad[36]={0};tinymekf_state before=state;CHECK(tinymekf_initialize(&state,q,zero,bad,&work)==TINYMEKF_INVALID_COVARIANCE);CHECK(!memcmp(&before,&state,sizeof state));}
 printf("tinymekf: 64 FilterPy predictions / 128 Joseph+reset updates; max quaternion %.3g, covariance %.3g\n",max_q,max_covariance);return 0;}
