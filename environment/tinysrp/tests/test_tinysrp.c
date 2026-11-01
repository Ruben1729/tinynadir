/* SPDX-License-Identifier: MIT */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL %s:%d %s\n",__FILE__,__LINE__,#c); return 1; } } while (0)

#include "tinysrp/tinysrp.h"
#include "data/reference.h"
int main(void){double out[3];
    for(size_t i=0;i<sizeof cases/sizeof cases[0];++i){const reference_case *c=&cases[i];CHECK(tinysrp_acceleration(c->r,c->sun,c->cr,c->area,c->mass,c->light,c->pressure,c->reference,out)==0);for(size_t j=0;j<3;++j)CHECK(fabs(out[j]-c->expected[j])<5e-15*fmax(1e-100,fabs(c->expected[j])));}
    {double r[3]={0},sun[3]={1,0,0};CHECK(tinysrp_acceleration(r,sun,1,1,1,1,1,1,out)==0&&out[0]==-1&&out[1]==0);CHECK(tinysrp_acceleration(r,sun,1,1,1,0,1,1,out)==0&&out[0]==0);out[0]=123;CHECK(tinysrp_acceleration(r,r,1,1,1,1,1,1,out)==TINYSRP_SINGULAR&&out[0]==123);CHECK(tinysrp_acceleration(r,sun,1,1,1,2,1,1,out)==TINYSRP_INVALID);}
    puts("tinysrp: independent inverse-square cannonball force cases pass");return 0;}
