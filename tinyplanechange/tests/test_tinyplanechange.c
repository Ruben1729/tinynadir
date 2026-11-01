/* SPDX-License-Identifier: MIT */
#include "tinyplanechange/tinyplanechange.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "data/reference.h"
#define CHECK(x) do { if (!(x)) {fprintf(stderr,"line %d: %s\n",__LINE__,#x);return 1;} } while(0)
int main(void)
{
    size_t i;double maximum=0;
    for(i=0;i<sizeof reference/sizeof reference[0];++i) {
        const double *r=reference[i];tinyplanechange_result out;
        CHECK(tinyplanechange_combined(r[0],r[1],r[2],r[3],r[4],&out)==TINYPLANECHANGE_OK);
        if(fabs(out.delta_v_m_s-r[5])>maximum)maximum=fabs(out.delta_v_m_s-r[5]);CHECK(fabs(out.delta_v_m_s-r[5])<1e-8);
        CHECK(out.unique_intersection);
        /* Either antipodal line-of-intersection direction is equivalent. */
        CHECK(fabs(remainder(out.initial_argument_latitude_rad-r[6],3.14159265358979323846))<1e-10);
        CHECK(fabs(remainder(out.final_argument_latitude_rad-r[7],3.14159265358979323846))<1e-10);
        if(r[2]==r[3]){tinyplanechange_result same;CHECK(tinyplanechange_node(r[0],r[1],r[2],r[4],&same)==TINYPLANECHANGE_OK);CHECK(fabs(same.delta_v_m_s-out.delta_v_m_s)<1e-12);}
    }
    for(i=0;i<sizeof incl_reference/sizeof incl_reference[0];++i){double dv;const double *r=incl_reference[i];CHECK(tinyplanechange_inclination(r[0],r[1],r[2],&dv)==TINYPLANECHANGE_OK);CHECK(fabs(dv-r[3])<1e-9);}
    {
        tinyplanechange_result out={1,2,3,4,5},saved=out;
        CHECK(tinyplanechange_combined(NAN,0,1,1,0,&out)==TINYPLANECHANGE_INVALID);CHECK(!memcmp(&out,&saved,sizeof out));
        CHECK(tinyplanechange_combined(7500,0,1,1,0,&out)==TINYPLANECHANGE_OK);CHECK(out.delta_v_m_s==0 && out.unique_intersection==0);
    }
    printf("Plane changes: 36 combined +27 inclination Vallado cases; max dv %.9g m/s\n",maximum);
    return 0;
}
