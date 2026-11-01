#include "tinycalibration/tinycalibration.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL %d: %s\n",__LINE__,#c);return 1;}}while(0)
#include "data/reference.h"
int main(void){double out[3],maximum=0;tinycalibration_parameters p;
 for(size_t c=0;c<sizeof reference_cases/sizeof reference_cases[0];++c){memcpy(p.bias,reference_cases[c].bias,sizeof p.bias);memcpy(p.correction,reference_cases[c].correction,sizeof p.correction);memcpy(p.sensor_to_body,reference_cases[c].alignment,sizeof p.sensor_to_body);CHECK(tinycalibration_apply(&p,reference_cases[c].raw,out)==TINYCALIBRATION_OK);for(size_t i=0;i<3;++i){double e=fabs(out[i]-reference_cases[c].expected[i]);CHECK(e<5e-14*fmax(1,fabs(reference_cases[c].expected[i])));maximum=fmax(maximum,e);}}
 {double raw[]={1,2,3},saved[]={123,456,789};memcpy(out,saved,sizeof out);p.sensor_to_body[0]=2;CHECK(tinycalibration_apply(&p,raw,out)==TINYCALIBRATION_INVALID_ALIGNMENT);CHECK(!memcmp(out,saved,sizeof out));CHECK(tinycalibration_apply(NULL,raw,out)==TINYCALIBRATION_INVALID_ARGUMENT);raw[0]=NAN;CHECK(tinycalibration_apply(&p,raw,out)==TINYCALIBRATION_INVALID_ARGUMENT);}
 printf("tinycalibration: 64 NumPy affine transforms; max error %.3g\n",maximum);return 0;}
