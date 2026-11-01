#include "tinyguidance/tinyguidance.h"
#include "tinyrotation/tinyrotation.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL %d: %s\n",__LINE__,#c);return 1;}}while(0)
#include "data/reference.h"
static double distance(const double a[4],const double b[4]){double m=0,p=0;for(size_t i=0;i<4;++i){m+=(a[i]-b[i])*(a[i]-b[i]);p+=(a[i]+b[i])*(a[i]+b[i]);}return sqrt(fmin(m,p));}
int main(void){double q[4],rate[3],maximum=0;
 for(size_t c=0;c<sizeof reference_cases/sizeof reference_cases[0];++c){CHECK(tinyguidance_nadir(reference_cases[c].position,reference_cases[c].velocity,q,rate)==TINYGUIDANCE_OK);double e=distance(q,reference_cases[c].expected);CHECK(e<3e-13);maximum=fmax(maximum,e);for(size_t i=0;i<3;++i)CHECK(fabs(rate[i]-reference_cases[c].rate[i])<3e-15);double nadir[3],body[3],norm=0;for(size_t i=0;i<3;++i)norm+=reference_cases[c].position[i]*reference_cases[c].position[i];for(size_t i=0;i<3;++i)nadir[i]=-reference_cases[c].position[i]/sqrt(norm);CHECK(tinyrotation_apply(q,nadir,body)==TINYROTATION_OK);CHECK(fabs(body[0])+fabs(body[1])+fabs(body[2]-1)<3e-13);}
 {double r[]={1,0,0},v[]={0,2,0};CHECK(tinyguidance_nadir(r,v,q,rate)==TINYGUIDANCE_OK);CHECK(rate[0]==0&&fabs(rate[1]+2)<1e-15&&rate[2]==0);}
 {double r[]={1,0,0},target[]={2,0,0},secondary[]={0,1,0},primary_body[]={0,0,1},secondary_body[]={1,0,0},mapped[3];CHECK(tinyguidance_target(r,target,secondary,primary_body,secondary_body,q)==TINYGUIDANCE_OK);CHECK(tinyrotation_apply(q,primary_body,mapped)==TINYROTATION_OK);double line[]={1,0,0};CHECK(tinyrotation_apply(q,line,mapped)==TINYROTATION_OK);CHECK(fabs(mapped[2]-1)<1e-14);}
 {double r[]={1,0,0},v[]={2,0,0},saved[]={123,456,789,42};memcpy(q,saved,sizeof q);CHECK(tinyguidance_nadir(r,v,q,rate)==TINYGUIDANCE_DEGENERATE);CHECK(!memcmp(q,saved,sizeof q));r[0]=NAN;CHECK(tinyguidance_nadir(r,v,q,rate)==TINYGUIDANCE_INVALID_ARGUMENT);}
 printf("tinyguidance: 64 independent NumPy/SciPy nadir frames; max quaternion error %.3g\n",maximum);return 0;}
