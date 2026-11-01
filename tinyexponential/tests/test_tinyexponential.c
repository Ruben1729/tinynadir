/* SPDX-License-Identifier: MIT */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL %s:%d %s\n",__FILE__,__LINE__,#c); return 1; } } while (0)

#include "tinyexponential/tinyexponential.h"
#include "data/reference.h"
int main(void){double out=123;
    for(size_t i=0;i<sizeof cases/sizeof cases[0];++i){const reference_case *c=&cases[i];CHECK(tinyexponential_density(c->layers,3,c->maximum,c->altitude,&out)==0);CHECK(fabs(out-c->expected)<2e-14*fmax(1e-300,fabs(c->expected)));}
    {tinyexponential_layer layer={0,1,1000};CHECK(tinyexponential_density(&layer,1,10000,0,&out)==0&&out==1);out=123;CHECK(tinyexponential_density(&layer,1,10000,-1,&out)==TINYEXPONENTIAL_RANGE&&out==123);layer.scale_height_m=0;CHECK(tinyexponential_density(&layer,1,10000,0,&out)==TINYEXPONENTIAL_INVALID);}
    puts("tinyexponential: independently evaluated layered density cases pass");return 0;}
