/* SPDX-License-Identifier: MIT */
#include "tinykepler/tinykepler.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "data/reference.h"
#define CHECK(x) do { if (!(x)) { fprintf(stderr,"line %d: %s\n",__LINE__,#x);return 1;} } while(0)
static double norm(const double x[3]) { return hypot(hypot(x[0],x[1]),x[2]); }
static double energy(double mu,const tinykepler_state *s)
{ return 0.5*(s->velocity_m_s[0]*s->velocity_m_s[0]+s->velocity_m_s[1]*s->velocity_m_s[1]+s->velocity_m_s[2]*s->velocity_m_s[2])-mu/norm(s->position_m); }
int main(void)
{
    size_t i,j;double max_r=0,max_v=0;
    for(i=0;i<sizeof reference/sizeof reference[0];++i) {
        tinykepler_state in,out,back;
        for(j=0;j<3;++j){in.position_m[j]=reference[i][2+j];in.velocity_m_s[j]=reference[i][5+j];}
        CHECK(tinykepler_propagate(reference[i][0],&in,reference[i][1],&out)==TINYKEPLER_OK);
        for(j=0;j<3;++j) {
            double er=fabs(out.position_m[j]-reference[i][8+j]),ev=fabs(out.velocity_m_s[j]-reference[i][11+j]);
            if(er>max_r)max_r=er;
            if(ev>max_v)max_v=ev;
            if(er>=1e-4 || ev>=1e-7) fprintf(stderr,"case %zu axis %zu: %.9g m %.9g m/s\n",i,j,er,ev);
            CHECK(er<1e-4);CHECK(ev<1e-7);
        }
        CHECK(fabs(energy(reference[i][0],&out)-energy(reference[i][0],&in))<2e-6);
        CHECK(tinykepler_propagate(reference[i][0],&out,-reference[i][1],&back)==TINYKEPLER_OK);
        for(j=0;j<3;++j){CHECK(fabs(back.position_m[j]-in.position_m[j])<1e-4);CHECK(fabs(back.velocity_m_s[j]-in.velocity_m_s[j])<1e-7);}
    }
    {
        tinykepler_state zero={{0,0,0},{0,0,0}},out={{1,2,3},{4,5,6}},saved=out;
        CHECK(tinykepler_propagate(1,&zero,1,&out)==TINYKEPLER_INVALID);
        CHECK(tinykepler_propagate(NAN,&out,1,&out)==TINYKEPLER_INVALID);
        CHECK(tinykepler_propagate(1,NULL,1,&out)==TINYKEPLER_INVALID);
        CHECK(!memcmp(&out,&saved,sizeof out));
        zero.position_m[0]=1;zero.velocity_m_s[0]=1;
        CHECK(tinykepler_propagate(1,&zero,1,&out)==TINYKEPLER_INVALID);
        CHECK(!memcmp(&out,&saved,sizeof out));
    }
    printf("Kepler: 42 DOP853 cases; max component error %.9g m / %.9g m/s\n",max_r,max_v);
    return 0;
}
