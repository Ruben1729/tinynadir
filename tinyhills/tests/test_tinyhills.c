/* SPDX-License-Identifier: MIT */
#include "tinyhills/tinyhills.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "data/reference.h"
#define CHECK(x) do { if (!(x)) {fprintf(stderr,"line %d: %s\n",__LINE__,#x);return 1;} } while(0)
int main(void)
{
    size_t i,j;double max_r=0,max_v=0,max_target=0;
    for(i=0;i<sizeof reference/sizeof reference[0];++i) {
        tinyhills_state in,out,back;double P[36];
        for(j=0;j<3;++j){in.position_m[j]=reference[i][2+j];in.velocity_m_s[j]=reference[i][5+j];}
        CHECK(tinyhills_propagate(reference[i][0],reference[i][1],&in,&out)==TINYHILLS_OK);
        CHECK(tinyhills_transition(reference[i][0],reference[i][1],P)==TINYHILLS_OK);
        for(j=0;j<3;++j){double er=fabs(out.position_m[j]-reference[i][8+j]),ev=fabs(out.velocity_m_s[j]-reference[i][11+j]);if(er>max_r)max_r=er;if(ev>max_v)max_v=ev;CHECK(er<2e-9);CHECK(ev<2e-12);}
        CHECK(tinyhills_propagate(reference[i][0],-reference[i][1],&out,&back)==TINYHILLS_OK);
        for(j=0;j<3;++j){CHECK(fabs(back.position_m[j]-in.position_m[j])<3e-8);CHECK(fabs(back.velocity_m_s[j]-in.velocity_m_s[j])<3e-11);}
    }
    for(i=0;i<sizeof target_reference/sizeof target_reference[0];++i) {
        tinyhills_state in,out;double v[3];
        CHECK(tinyhills_rendezvous(target_reference[i][0],target_reference[i][1],target_reference[i]+2,v)==TINYHILLS_OK);
        for(j=0;j<3;++j){double error=fabs(v[j]-target_reference[i][5+j]);if(error>max_target)max_target=error;CHECK(error<2e-11);in.position_m[j]=target_reference[i][2+j];in.velocity_m_s[j]=v[j];}
        CHECK(tinyhills_propagate(target_reference[i][0],target_reference[i][1],&in,&out)==TINYHILLS_OK);
        for(j=0;j<3;++j)CHECK(fabs(out.position_m[j])<1e-8);
    }
    {
        double r[3]={1,2,3},v[3]={7,8,9},saved[3]={7,8,9};
        CHECK(tinyhills_rendezvous(1,3.14159265358979323846,r,v)==TINYHILLS_SINGULAR);CHECK(!memcmp(v,saved,sizeof v));
        CHECK(tinyhills_rendezvous(0,2,r,v)==TINYHILLS_OK);CHECK(v[0]==-.5 && v[1]==-1 && v[2]==-1.5);
        CHECK(tinyhills_rendezvous(NAN,2,r,v)==TINYHILLS_INVALID);
    }
    printf("Hills: 30 propagation +24 targeting cases; max %.9g m / %.9g m/s; targeting %.9g m/s\n",max_r,max_v,max_target);
    return 0;
}
