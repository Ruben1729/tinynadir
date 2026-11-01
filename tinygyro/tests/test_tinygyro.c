#include "tinygyro/tinygyro.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL %d: %s\n",__LINE__,#c);return 1;}}while(0)
#include "data/reference.h"
static double distance(const double a[4],const double b[4]){double m=0,p=0;for(size_t i=0;i<4;++i){m+=(a[i]-b[i])*(a[i]-b[i]);p+=(a[i]+b[i])*(a[i]+b[i]);}return sqrt(fmin(m,p));}
int main(void){double out[4],derivative[4],maximum=0;
 for(size_t c=0;c<sizeof reference_cases/sizeof reference_cases[0];++c){CHECK(tinygyro_step(reference_cases[c].q,reference_cases[c].rate,reference_cases[c].bias,reference_cases[c].dt,out)==TINYGYRO_OK);double e=distance(out,reference_cases[c].expected);CHECK(e<3e-14);maximum=fmax(maximum,e);CHECK(tinygyro_derivative(reference_cases[c].q,reference_cases[c].rate,reference_cases[c].bias,derivative)==TINYGYRO_OK);for(size_t i=0;i<4;++i)CHECK(fabs(derivative[i]-reference_cases[c].derivative[i])<3e-10);}
 {double q[]={1,0,0,0},rate[]={0,0,1};CHECK(tinygyro_step(q,rate,NULL,.2,out)==TINYGYRO_OK);CHECK(fabs(out[0]-cos(.1))<1e-15&&fabs(out[3]+sin(.1))<1e-15);CHECK(tinygyro_derivative(q,rate,NULL,derivative)==TINYGYRO_OK);CHECK(derivative[3]==-.5);}
 {double q[]={1,0,0,0},rate[]={0,0,1},saved[4]={123,456,789,42};memcpy(out,saved,sizeof out);CHECK(tinygyro_step(q,rate,NULL,-1,out)==TINYGYRO_INVALID_ARGUMENT);CHECK(!memcmp(out,saved,sizeof out));rate[0]=NAN;CHECK(tinygyro_step(q,rate,NULL,1,out)==TINYGYRO_INVALID_ARGUMENT);CHECK(tinygyro_step(NULL,rate,NULL,1,out)==TINYGYRO_INVALID_ARGUMENT);}
 printf("tinygyro: 64 SciPy held-rate rotations; max quaternion error %.3g\n",maximum);return 0;}
