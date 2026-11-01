/* SPDX-License-Identifier: MIT */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL %s:%d %s\n",__FILE__,__LINE__,#c); return 1; } } while (0)

#include "tinydrag/tinydrag.h"
#include "data/reference.h"
int main(void){double out[3];
    for(size_t i=0;i<sizeof cases/sizeof cases[0];++i){const reference_case *c=&cases[i];CHECK(tinydrag_acceleration(c->density,c->velocity,c->cd,c->area,c->mass,out)==0);for(size_t j=0;j<3;++j)CHECK(fabs(out[j]-c->expected[j])<5e-15*fmax(1e-100,fabs(c->expected[j])));CHECK(out[0]*c->velocity[0]+out[1]*c->velocity[1]+out[2]*c->velocity[2]<=0);}
    {double v[3]={1,2,3};out[0]=123;CHECK(tinydrag_acceleration(-1,v,2,1,1,out)==TINYDRAG_INVALID&&out[0]==123);CHECK(tinydrag_acceleration(0,v,2,1,1,out)==0&&out[0]==0);CHECK(tinydrag_acceleration(1,v,2,1,0,out)==TINYDRAG_INVALID);v[1]=NAN;CHECK(tinydrag_acceleration(1,v,2,1,1,out)==TINYDRAG_INVALID);}
    puts("tinydrag: independent vector drag and dissipative-work checks pass");return 0;}
