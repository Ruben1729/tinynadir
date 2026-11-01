#include "tinyslew/tinyslew.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL %d: %s\n",__LINE__,#c);return 1;}}while(0)
#include "data/reference.h"
static double distance(const double a[4],const double b[4]){double m=0,p=0;for(size_t i=0;i<4;++i){m+=(a[i]-b[i])*(a[i]-b[i]);p+=(a[i]+b[i])*(a[i]+b[i]);}return sqrt(fmin(m,p));}
int main(void){double q[4],rate[3],accel[3],maximum=0;tinyslew_plan plan;
 for(size_t c=0;c<sizeof reference_cases/sizeof reference_cases[0];++c){CHECK(tinyslew_initialize(&plan,reference_cases[c].initial,reference_cases[c].target,reference_cases[c].vmax,reference_cases[c].amax)==TINYSLEW_OK);CHECK(fabs(plan.duration-reference_cases[c].duration)<3e-13*fmax(1,plan.duration));for(size_t s=0;s<8;++s){double t=reference_cases[c].sample[s].t;CHECK(tinyslew_sample(&plan,t,q,rate,accel)==TINYSLEW_OK);double e=distance(q,reference_cases[c].sample[s].q);CHECK(e<3e-13);maximum=fmax(maximum,e);double norm=0;for(size_t i=0;i<3;++i){CHECK(fabs(rate[i]-reference_cases[c].sample[s].rate[i])<3e-13);norm+=rate[i]*rate[i];if(fabs(t-plan.accel_time)>1e-12&&fabs(t-plan.accel_time-plan.cruise_time)>1e-12&&fabs(t-plan.duration)>1e-12)CHECK(fabs(accel[i]-reference_cases[c].sample[s].acceleration[i])<3e-13);}CHECK(sqrt(norm)<=reference_cases[c].vmax+1e-13);}}
 {double initial[]={1,0,0,0},target[]={-1,0,0,0};CHECK(tinyslew_initialize(&plan,initial,target,.1,.01)==TINYSLEW_OK);CHECK(plan.duration==0);CHECK(tinyslew_sample(&plan,0,q,rate,accel)==TINYSLEW_OK);CHECK(q[0]==1&&rate[0]==0&&accel[0]==0);}
 {double initial[]={1,0,0,0},target[]={cos(.1),0,0,-sin(.1)};CHECK(tinyslew_initialize(&plan,initial,target,.1,.01)==TINYSLEW_OK);CHECK(tinyslew_sample(&plan,plan.duration*.25,q,rate,accel)==TINYSLEW_OK);CHECK(rate[2]>0);}
 {double initial[]={1,0,0,0};tinyslew_plan saved=plan;CHECK(tinyslew_initialize(&plan,initial,initial,0,.1)==TINYSLEW_INVALID_ARGUMENT);CHECK(!memcmp(&plan,&saved,sizeof plan));double savedq[]={123,456,789,42};memcpy(q,savedq,sizeof q);CHECK(tinyslew_sample(&plan,-1,q,rate,accel)==TINYSLEW_INVALID_ARGUMENT);CHECK(!memcmp(q,savedq,sizeof q));}
 printf("tinyslew: 32 independent profiles / 256 SciPy SLERP samples; max quaternion error %.3g\n",maximum);return 0;}
