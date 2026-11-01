/* SPDX-License-Identifier: MIT */
#include "tinyanomaly/tinyanomaly.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "data/reference.h"
#define CHECK(x) do { if (!(x)) { fprintf(stderr,"line %d: %s\n",__LINE__,#x);return 1;} } while(0)
int main(void)
{
    size_t i;double max_a=0,max_nu=0;
    for(i=0;i<sizeof reference/sizeof reference[0];++i) {
        tinyanomaly_solution out;
        CHECK(tinyanomaly_from_mean(reference[i][0],reference[i][1],&out)==TINYANOMALY_OK);
        if(fabs(out.anomaly_rad-reference[i][2])>max_a)max_a=fabs(out.anomaly_rad-reference[i][2]);
        if(fabs(out.true_anomaly_rad-reference[i][3])>max_nu)max_nu=fabs(out.true_anomaly_rad-reference[i][3]);
        CHECK(fabs(out.anomaly_rad-reference[i][2])<2e-12);
        CHECK(fabs(out.true_anomaly_rad-reference[i][3])<2e-10);
    }
    for(i=0;i<sizeof inverse_reference/sizeof inverse_reference[0];++i) {
        double m,a;
        CHECK(tinyanomaly_from_true(inverse_reference[i][0],inverse_reference[i][1],&m,&a)==TINYANOMALY_OK);
        CHECK(fabs(m-inverse_reference[i][2])<2e-11*fmax(1.,fabs(m)));
        CHECK(fabs(a-inverse_reference[i][3])<2e-12*fmax(1.,fabs(a)));
    }
    {
        tinyanomaly_solution out={12,13},saved=out;double m=9,a=8;
        CHECK(tinyanomaly_from_mean(-1,1,&out)==TINYANOMALY_INVALID);
        CHECK(tinyanomaly_from_mean(1,NAN,&out)==TINYANOMALY_INVALID);
        CHECK(tinyanomaly_from_mean(0,0,NULL)==TINYANOMALY_INVALID);
        CHECK(!memcmp(&out,&saved,sizeof out));
        CHECK(tinyanomaly_from_true(2,3,&m,&a)==TINYANOMALY_INVALID);
        CHECK(m==9 && a==8);
    }
    printf("Anomaly: 80 mean and 35 inverse cases; max errors %.9g / %.9g rad\n",max_a,max_nu);
    return 0;
}
