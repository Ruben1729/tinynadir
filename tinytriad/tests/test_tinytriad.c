/* SPDX-License-Identifier: MIT */
#include "tinytriad/tinytriad.h"
#include "tinyrotation/tinyrotation.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL %s:%d: %s\n",__FILE__,__LINE__,#c);return 1;}}while(0)
#include "data/reference.h"
int main(void){
    double matrix[9],q[4],max_error=0;
    for(size_t c=0;c<sizeof reference_cases/sizeof reference_cases[0];++c){
        CHECK(tinytriad_dcm(reference_cases[c].r1,reference_cases[c].r2,reference_cases[c].b1,reference_cases[c].b2,matrix)==TINYTRIAD_OK);
        for(size_t i=0;i<9;++i){double e=fabs(matrix[i]-reference_cases[c].matrix[i]);CHECK(e<3e-14);max_error=fmax(max_error,e);}
        CHECK(tinytriad_solve(reference_cases[c].r1,reference_cases[c].r2,reference_cases[c].b1,reference_cases[c].b2,q)==TINYTRIAD_OK);
        CHECK(tinyrotation_to_dcm(q,matrix)==TINYROTATION_OK);
        for(size_t i=0;i<9;++i)CHECK(fabs(matrix[i]-reference_cases[c].matrix[i])<3e-14);
        for(size_t i=0;i<3;++i){double v=0;for(size_t j=0;j<3;++j)v+=matrix[i*3+j]*reference_cases[c].r1[j];CHECK(fabs(v-reference_cases[c].b1[i])<3e-14);}
    }
    {double a[3]={1,0,0},b[3]={0,1,0},zero[3]={0},near[3]={1,1e-12,0},bad[3]={NAN,1,0},out[9],saved[9];
      for(size_t i=0;i<9;++i)out[i]=saved[i]=123+i;
      CHECK(tinytriad_dcm(a,a,a,b,out)==TINYTRIAD_DEGENERATE);CHECK(tinytriad_dcm(a,near,a,b,out)==TINYTRIAD_DEGENERATE);
      CHECK(tinytriad_dcm(a,b,zero,b,out)==TINYTRIAD_DEGENERATE);CHECK(tinytriad_dcm(a,b,a,a,out)==TINYTRIAD_DEGENERATE);
      CHECK(tinytriad_dcm(bad,b,a,b,out)==TINYTRIAD_INVALID_ARGUMENT);CHECK(tinytriad_dcm(NULL,b,a,b,out)==TINYTRIAD_INVALID_ARGUMENT);
      CHECK(tinytriad_dcm(a,b,a,b,NULL)==TINYTRIAD_INVALID_ARGUMENT);CHECK(!memcmp(out,saved,sizeof out));}
    {double a[3]={1e300,0,0},b[3]={0,1e-300,0};CHECK(tinytriad_dcm(a,b,b,a,matrix)==TINYTRIAD_OK);CHECK(fabs(matrix[1]-1)<1e-15);}
    printf("tinytriad: 128 AHRS cases; max DCM difference %.3g\n",max_error);return 0;
}
