/* SPDX-License-Identifier: MIT */
#include "tinyelements/tinyelements.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "data/reference.h"
#define CHECK(x) do { if (!(x)) { fprintf(stderr,"line %d: %s\n",__LINE__,#x);return 1;} } while(0)
int main(void)
{
    size_t i,j;double max_r=0,max_v=0,max_angle=0;
    for(i=0;i<sizeof reference/sizeof reference[0];++i) {
        const double *ref=reference[i];
        tinyelements_classical e={ref[1],ref[2],ref[3],ref[4],ref[5],ref[6]},back;
        tinyelements_state s,restored;tinyelements_equinoctial q;double C[9];
        CHECK(tinyelements_from_classical(ref[0],&e,&s)==TINYELEMENTS_OK);
        for(j=0;j<3;++j) {
            double er=fabs(s.position_m[j]-ref[7+j]),ev=fabs(s.velocity_m_s[j]-ref[10+j]);
            if(er>max_r)max_r=er;
            if(ev>max_v)max_v=ev;
            CHECK(er<2e-7);CHECK(ev<2e-10);
        }
        CHECK(tinyelements_to_classical(ref[0],&s,&back)==TINYELEMENTS_OK);
        CHECK(fabs(back.p_m-ref[13])<2e-7);CHECK(fabs(back.eccentricity-ref[14])<1e-13);
        {double a[]={back.inclination_rad,back.raan_rad,back.arg_periapsis_rad,back.true_anomaly_rad};for(j=0;j<4;++j){double error=fabs(remainder(a[j]-ref[15+j],6.2831853071795864769));if(error>max_angle)max_angle=error;CHECK(error<1e-12);}}
        CHECK(tinyelements_to_equinoctial(ref[0],&s,&q)==TINYELEMENTS_OK);
        CHECK(tinyelements_from_equinoctial(ref[0],&q,&restored)==TINYELEMENTS_OK);
        for(j=0;j<3;++j){CHECK(fabs(restored.position_m[j]-s.position_m[j])<3e-7);CHECK(fabs(restored.velocity_m_s[j]-s.velocity_m_s[j])<3e-10);}
        CHECK(tinyelements_local_frame(ref[0],&s,TINYELEMENTS_RSW,C)==TINYELEMENTS_OK);
        for(j=0;j<9;++j)CHECK(fabs(C[j]-ref[19+j])<1e-14);
        CHECK(tinyelements_local_frame(ref[0],&s,TINYELEMENTS_NTW,C)==TINYELEMENTS_OK);
        for(j=0;j<9;++j)CHECK(fabs(C[j]-ref[28+j])<1e-14);
        CHECK(tinyelements_local_frame(ref[0],&s,TINYELEMENTS_PQW,C)==TINYELEMENTS_OK);
        for(j=0;j<3;++j){size_t k;double r=0,v=0;for(k=0;k<3;++k){r+=C[3*j+k]*s.position_m[k];v+=C[3*j+k]*s.velocity_m_s[k];}CHECK(fabs(r-ref[37+j])<3e-7);CHECK(fabs(v-ref[40+j])<3e-10);}
    }
    for(i=0;i<6;++i) {
        tinyelements_classical e={7000000,i%2 ? .1 : 0,i<2 ? 0 : i<4 ? 1.1 : 3.14159265358979323846,.8,1.1,.3},back;
        tinyelements_state a,b;
        CHECK(tinyelements_from_classical(398600.4415e9,&e,&a)==TINYELEMENTS_OK);
        CHECK(tinyelements_to_classical(398600.4415e9,&a,&back)==TINYELEMENTS_OK);
        CHECK(tinyelements_from_classical(398600.4415e9,&back,&b)==TINYELEMENTS_OK);
        for(j=0;j<3;++j){CHECK(fabs(a.position_m[j]-b.position_m[j])<2e-7);CHECK(fabs(a.velocity_m_s[j]-b.velocity_m_s[j])<2e-10);}
    }
    {
        tinyelements_state a={{1,0,0},{1,0,0}};tinyelements_classical e={1,2,3,4,5,6},saved=e;
        CHECK(tinyelements_to_classical(1,&a,&e)==TINYELEMENTS_SINGULAR);CHECK(!memcmp(&e,&saved,sizeof e));
        CHECK(tinyelements_to_classical(NAN,&a,&e)==TINYELEMENTS_INVALID);
    }
    printf("Elements: 15 Vallado states and local bases; max %.9g m / %.9g m/s, angles %.9g rad\n",max_r,max_v,max_angle);
    return 0;
}
