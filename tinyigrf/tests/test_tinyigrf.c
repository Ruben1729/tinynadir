#include "tinyigrf/tinyigrf.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(x) do{if(!(x)){fprintf(stderr,"%s:%d %s\n",__FILE__,__LINE__,#x);exit(1);}}while(0)
#include "data/reference.h"
#include "data/geodetic.h"
int main(void) {
    size_t i,k;double max_error=0,b[3],pole[3],near[3],r[3],original[3]={1,2,3};
    for(i=0;i<sizeof references/sizeof references[0];i++) {
        CHECK(tinyigrf_geocentric(references[i].year,references[i].radius,references[i].theta,references[i].phi,references[i].degree,b)==TINYIGRF_OK);
        for(k=0;k<3;k++){double e=fabs(b[k]-references[i].b[k]);CHECK(e<5e-14);if(e>max_error)max_error=e;}
    }
    for(i=0;i<sizeof geodetic_cases/sizeof geodetic_cases[0];i++) {
        CHECK(tinyigrf_geodetic(geodetic_cases[i][0],geodetic_cases[i][1],geodetic_cases[i][2],geodetic_cases[i][3],13,b)==TINYIGRF_OK);
        for(k=0;k<3;k++)CHECK(fabs(b[k]-geodetic_cases[i][4+k])<1e-12);
    }
    /* Exact north-pole spherical limit, where the reference package divides by zero. */
    CHECK(tinyigrf_geocentric(2027,6771200,0,.7,13,pole)==TINYIGRF_OK);
    CHECK(tinyigrf_geocentric(2027,6771200,1e-11,.7,13,near)==TINYIGRF_OK);
    for(k=0;k<3;k++)CHECK(fabs(pole[k]-near[k])<1e-14);
    r[0]=0;r[1]=0;r[2]=6771200;
    CHECK(tinyigrf_ecef(2027,r,13,b)==TINYIGRF_OK);CHECK(isfinite(b[0])&&isfinite(b[1])&&isfinite(b[2]));
    CHECK(tinyigrf_ecef(2027,r,13,r)==TINYIGRF_OK);for(k=0;k<3;k++)CHECK(r[k]==b[k]);
    memcpy(b,original,sizeof b);
    CHECK(tinyigrf_geocentric(2024,6771200,1,.7,13,b)==TINYIGRF_RANGE);CHECK(!memcmp(b,original,sizeof b));
    CHECK(tinyigrf_geocentric(2030.001,6771200,1,.7,13,b)==TINYIGRF_RANGE);
    CHECK(tinyigrf_geocentric(2027,6771200,1,.7,14,b)==TINYIGRF_RANGE);
    CHECK(tinyigrf_geocentric(2027,6771200,1,.7,0,b)==TINYIGRF_RANGE);
    CHECK(tinyigrf_geocentric(NAN,6771200,1,.7,13,b)==TINYIGRF_INVALID);
    CHECK(tinyigrf_geocentric(2027,0,1,.7,13,b)==TINYIGRF_RANGE);
    CHECK(tinyigrf_geocentric(2027,6771200,1,.7,13,NULL)==TINYIGRF_INVALID);
    CHECK(!memcmp(b,original,sizeof b));
    printf("tinyigrf: %zu ppigrf cases, geodetic/ECEF and pole checks pass; max component error %.9g nT\n",sizeof references/sizeof references[0],max_error*1e9);
    return 0;
}
