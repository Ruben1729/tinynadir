#include "tinycomplementary/tinycomplementary.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL %d: %s\n",__LINE__,#c);return 1;}}while(0)
#include "data/reference.h"
static double distance(const double a[4],const double b[4]){double m=0,p=0;for(size_t i=0;i<4;++i){m+=(a[i]-b[i])*(a[i]-b[i]);p+=(a[i]+b[i])*(a[i]+b[i]);}return sqrt(fmin(m,p));}
int main(void){tinycomplementary_state state;tinycomplementary_parameters p={2,.2,.1};const double zero[3]={0},reference[3]={0,0,1};double maximum=0;
 CHECK(tinycomplementary_initialize(&state,initial_quaternion,zero)==TINYCOMPLEMENTARY_OK);
 for(size_t c=0;c<sizeof reference_cases/sizeof reference_cases[0];++c){CHECK(tinycomplementary_step(&state,&p,reference_cases[c].gyro,.01,1,reference,reference_cases[c].body,NULL)==TINYCOMPLEMENTARY_OK);double e=distance(state.q,reference_cases[c].expected);CHECK(e<3e-13);maximum=fmax(maximum,e);for(size_t i=0;i<3;++i)CHECK(fabs(state.bias[i]-reference_cases[c].bias[i])<3e-13);}
 {tinycomplementary_state before=state;CHECK(tinycomplementary_step(&state,&p,zero,0,1,reference,reference,NULL)==TINYCOMPLEMENTARY_INVALID_ARGUMENT);CHECK(!memcmp(&state,&before,sizeof state));CHECK(tinycomplementary_step(&state,&p,zero,.01,65,reference,reference,NULL)==TINYCOMPLEMENTARY_INVALID_ARGUMENT);}
 {double q[]={cos(.2),0,0,sin(.2)};tinycomplementary_state test;CHECK(tinycomplementary_initialize(&test,q,zero)==TINYCOMPLEMENTARY_OK);for(size_t i=0;i<100;++i)CHECK(tinycomplementary_step(&test,&p,zero,.01,1,reference,reference,NULL)==TINYCOMPLEMENTARY_OK);CHECK(distance(test.q,q)<1e-14);}
 {double q[]={cos(.2),sin(.2),0,0},refs[]={1,0,0,0,0,1};tinycomplementary_state test;tinycomplementary_parameters cfg={2,0,.1};CHECK(tinycomplementary_initialize(&test,q,zero)==TINYCOMPLEMENTARY_OK);for(size_t i=0;i<1000;++i)CHECK(tinycomplementary_step(&test,&cfg,zero,.01,2,refs,refs,NULL)==TINYCOMPLEMENTARY_OK);CHECK(fabs(test.q[1])<1e-5);}
 {double q[]={1,0,0,0},bias[]={0,0,.5};tinycomplementary_state test;CHECK(tinycomplementary_initialize(&test,q,bias)==TINYCOMPLEMENTARY_OK);CHECK(tinycomplementary_step(&test,&p,zero,.01,0,NULL,NULL,NULL)==TINYCOMPLEMENTARY_OK);CHECK(test.bias[2]==.1);}
 printf("tinycomplementary: 128 AHRS Mahony samples; max quaternion difference %.3g\n",maximum);return 0;}
