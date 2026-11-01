/* SPDX-License-Identifier: MIT */
#include "tinyhohmann/tinyhohmann.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "data/reference.h"
#define CHECK(x) do { if (!(x)) {fprintf(stderr,"line %d: %s\n",__LINE__,#x);return 1;} } while(0)
int main(void)
{
    size_t i,j;double max_dv=0,max_time=0;
    for(i=0;i<sizeof reference/sizeof reference[0];++i){tinyhohmann_transfer t;CHECK(tinyhohmann_circular(reference[i][0],reference[i][1],reference[i][2],&t)==TINYHOHMANN_OK);for(j=0;j<2;++j){double er=fabs(t.burn_m_s[j]-reference[i][3+j]);if(er>max_dv)max_dv=er;CHECK(er<1e-8);}CHECK(t.burn_count==2);CHECK(fabs(t.time_of_flight_s-reference[i][5])<1e-7);}
    for(i=0;i<sizeof bi_reference/sizeof bi_reference[0];++i){tinyhohmann_transfer t;CHECK(tinyhohmann_bielliptic(bi_reference[i][0],bi_reference[i][1],bi_reference[i][2],bi_reference[i][3],&t)==TINYHOHMANN_OK);for(j=0;j<3;++j){double er=fabs(t.burn_m_s[j]-bi_reference[i][4+j]);if(er>max_dv)max_dv=er;CHECK(er<1e-8);}CHECK(t.burn_count==3);CHECK(fabs(t.time_of_flight_s-bi_reference[i][7])<1e-6);}
    for(i=0;i<sizeof tangent_reference/sizeof tangent_reference[0];++i){tinyhohmann_tangent t;const double *r=tangent_reference[i];CHECK(tinyhohmann_one_tangent(r[0],r[1],r[2],(int)r[3],r[4],&t)==TINYHOHMANN_OK);for(j=0;j<2;++j){double er=fabs(t.transfer.burn_m_s[j]-r[5+j]);if(er>max_dv)max_dv=er;CHECK(er<1e-7);}if(fabs(t.transfer.time_of_flight_s-r[7])>max_time)max_time=fabs(t.transfer.time_of_flight_s-r[7]);CHECK(fabs(t.transfer.time_of_flight_s-r[7])<1e-6);CHECK(fabs(t.eccentricity-r[8])<1e-13);CHECK(fabs(t.p_m-r[9])<1e-6);}
    {
        tinyhohmann_transfer t={{1,2,3},4,5,6},saved=t;
        CHECK(tinyhohmann_circular(1,-1,2,&t)==TINYHOHMANN_INVALID);CHECK(!memcmp(&t,&saved,sizeof t));
        CHECK(tinyhohmann_bielliptic(1,2,1,3,&t)==TINYHOHMANN_GEOMETRY);CHECK(!memcmp(&t,&saved,sizeof t));
    }
    printf("Transfers: 5 Hohmann +15 bielliptic +6 one-tangent; max dv %.9g m/s; tangent time %.9g s\n",max_dv,max_time);
    return 0;
}
