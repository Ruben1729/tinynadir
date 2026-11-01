#include "tinydesaturation/tinydesaturation.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL %d: %s\n",__LINE__,#c);return 1;}}while(0)
#include "data/reference.h"
int main(void){double dipole[3],torque[3],residual[3],maximum=0;
 for(size_t c=0;c<sizeof reference_cases/sizeof reference_cases[0];++c){CHECK(tinydesaturation_command(reference_cases[c].momentum,reference_cases[c].target,reference_cases[c].gain,reference_cases[c].field,reference_cases[c].limit,dipole,torque,residual)==TINYDESATURATION_OK);double decay=0;for(size_t i=0;i<3;++i){double e=fabs(dipole[i]-reference_cases[c].expected[i]);CHECK(e<3e-13);CHECK(fabs(torque[i]-reference_cases[c].torque[i])<3e-17);CHECK(fabs(residual[i]-reference_cases[c].residual[i])<3e-17);maximum=fmax(maximum,e);decay+=(reference_cases[c].momentum[i]-reference_cases[c].target[i])*torque[i];}CHECK(decay<=1e-18);}
 {double momentum[]={0,0,.01},field[]={1e-5,0,0},limit[]={10,10,10};CHECK(tinydesaturation_command(momentum,NULL,.001,field,limit,dipole,torque,NULL)==TINYDESATURATION_OK);CHECK(torque[2]<0);CHECK(torque[0]==0);}
 {double momentum[]={1,0,0},field[]={0,0,0},limit[]={1,1,1},saved[]={123,456,789};memcpy(dipole,saved,sizeof dipole);memcpy(torque,saved,sizeof torque);CHECK(tinydesaturation_command(momentum,NULL,.1,field,limit,dipole,torque,NULL)==TINYDESATURATION_DEGENERATE);CHECK(!memcmp(dipole,saved,sizeof dipole)&&!memcmp(torque,saved,sizeof torque));}
 printf("tinydesaturation: 64 NumPy SVD unloading cases; max dipole error %.3g A m^2\n",maximum);return 0;}
