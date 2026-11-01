#include "tinymagnetic/tinymagnetic.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL %d: %s\n",__LINE__,#c);return 1;}}while(0)
#include "data/reference.h"
int main(void){double dipole[3],residual[3],maximum=0;
 for(size_t c=0;c<sizeof reference_cases/sizeof reference_cases[0];++c){CHECK(tinymagnetic_allocate(reference_cases[c].field,reference_cases[c].torque,reference_cases[c].limit,dipole,residual)==TINYMAGNETIC_OK);for(size_t i=0;i<3;++i){double e=fabs(dipole[i]-reference_cases[c].expected[i]);CHECK(e<3e-13);CHECK(fabs(residual[i]-reference_cases[c].residual[i])<3e-17);maximum=fmax(maximum,e);CHECK(fabs(dipole[i])<=reference_cases[c].limit[i]);}}
 {double field[]={0,0,2e-5},torque[]={0,0,1e-5},limit[]={1,1,1};CHECK(tinymagnetic_allocate(field,torque,limit,dipole,residual)==TINYMAGNETIC_OK);CHECK(dipole[0]==0&&dipole[1]==0&&dipole[2]==0&&residual[2]==torque[2]);}
 {double zero[]={0,0,0},limit[]={1,1,1},saved[]={123,456,789};memcpy(dipole,saved,sizeof dipole);memcpy(residual,saved,sizeof residual);CHECK(tinymagnetic_allocate(zero,zero,limit,dipole,residual)==TINYMAGNETIC_DEGENERATE);CHECK(!memcmp(dipole,saved,sizeof dipole)&&!memcmp(residual,saved,sizeof residual));limit[0]=-1;CHECK(tinymagnetic_allocate(zero,zero,limit,dipole,residual)==TINYMAGNETIC_INVALID_ARGUMENT);}
 printf("tinymagnetic: 64 NumPy SVD pseudoinverse cases; max dipole error %.3g A m^2\n",maximum);return 0;}
