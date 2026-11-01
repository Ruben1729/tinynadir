/* SPDX-License-Identifier: MIT */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL %s:%d %s\n",__FILE__,__LINE__,#c); return 1; } } while (0)

#include "tinythirdbody/tinythirdbody.h"
#include "data/reference.h"
int main(void){double out[3];double maximum=0;
    for(size_t i=0;i<sizeof cases/sizeof cases[0];++i){const reference_case *c=&cases[i];CHECK(tinythirdbody_acceleration(c->r,c->b,c->mu,out)==0);for(size_t j=0;j<3;++j){double e=fabs(out[j]-c->expected[j]);maximum=fmax(maximum,e);CHECK(e<2e-13*fmax(1e-12,fabs(c->expected[j])));}}
    {double r[3]={0,0,0},b[3]={1,0,0};CHECK(tinythirdbody_acceleration(r,b,1,out)==0&&out[0]==0);out[0]=123;CHECK(tinythirdbody_acceleration(b,b,1,out)==TINYTHIRDBODY_SINGULAR&&out[0]==123);CHECK(tinythirdbody_acceleration(r,b,-1,out)==TINYTHIRDBODY_INVALID);r[0]=NAN;CHECK(tinythirdbody_acceleration(r,b,1,out)==TINYTHIRDBODY_INVALID);}
    printf("tinythirdbody: high-precision differential gravity max absolute error %.12g m/s2\n",maximum);return 0;}
