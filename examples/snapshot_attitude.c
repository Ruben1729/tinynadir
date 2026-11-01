/* SPDX-License-Identifier: MIT */
#include "tinyomm/tinyomm.h"
#include "tinytime/tinytime.h"
#include "tinyframes/tinyframes.h"
#include "tinysun/tinysun.h"
#include "tinyigrf/tinyigrf.h"
#include "tinyrotation/tinyrotation.h"
#include "tinytriad/tinytriad.h"
#include <math.h>
#include <stdio.h>
#define REQUIRE(c) do { if (!(c)) { fprintf(stderr,"Example failed at line %d\n",__LINE__); return 1; } } while(0)
int main(void) {
    /* Illustrative elements and synthetic sensors, not an operational TLE. */
    const char omm[] =
        "CENTER_NAME = EARTH\nREF_FRAME = TEME\nTIME_SYSTEM = UTC\n"
        "MEAN_ELEMENT_THEORY = SGP4\nEPOCH = 2026-11-01T12:00:00Z\n"
        "MEAN_MOTION = 15\nECCENTRICITY = .001\nINCLINATION = 51\n"
        "RA_OF_ASC_NODE = 40\nARG_OF_PERICENTER = 20\nMEAN_ANOMALY = 30\nBSTAR = 0\n";
    tinytle_elements elements;
    tgp4_satellite satellite = {0};
    tgp4_state propagated;
    tinytime_jd utc,tai,tt,ut1,tdb;
    tinyframes_epoch epoch = {{0,0},{0,0},0,0,0,0,0};
    tinysun_state sun;
    double rt[3],vt[3],rg[3],vg[3],re[3],ve[3],magnetic_ecef[3],magnetic_ref[3];
    double matrix[3][3],sun_ref[3],magnetic_body[3],sun_body[3],truth[4],attitude[4],year;
    const double axis[3] = {1,2,3};
    REQUIRE(tinyomm_parse_kvn(omm,sizeof omm-1,0,&elements)==TINYOMM_OK);
    REQUIRE(tinytle_initialize_elements(&elements,&satellite)==TINYTLE_OK);
    REQUIRE(tgp4_propagate(&satellite,0,&propagated)==TGP4_OK);
    utc.first=elements.epoch_jd1;utc.fraction=elements.epoch_jd2;
    /* Example time/EOP offsets; flight applications supply their valid values. */
    REQUIRE(tinytime_utc_to_tai(utc,37,&tai)==TINYTIME_OK);
    REQUIRE(tinytime_tai_to_tt(tai,&tt)==TINYTIME_OK);
    REQUIRE(tinytime_utc_to_ut1(utc,0,&ut1)==TINYTIME_OK);
    REQUIRE(tinytime_tt_to_tdb(tt,&tdb)==TINYTIME_OK);
    REQUIRE(tinytime_decimal_year(utc,&year)==TINYTIME_OK);
    epoch.tt[0]=tt.first;epoch.tt[1]=tt.fraction;epoch.ut1[0]=ut1.first;epoch.ut1[1]=ut1.fraction;
    for(size_t i=0;i<3;++i){rt[i]=1000*propagated.position_km[i];vt[i]=1000*propagated.velocity_km_s[i];}
    REQUIRE(tinyframes_transform(&epoch,TINYFRAMES_TEME,TINYFRAMES_GCRF,TINYFRAMES_IAU2006A,rt,vt,rg,vg)==TINYFRAMES_OK);
    REQUIRE(tinyframes_transform(&epoch,TINYFRAMES_TEME,TINYFRAMES_ITRF,TINYFRAMES_IAU2006A,rt,vt,re,ve)==TINYFRAMES_OK);
    REQUIRE(tinyigrf_ecef(year,re,13,magnetic_ecef)==TINYIGRF_OK);
    REQUIRE(tinyframes_matrix(&epoch,TINYFRAMES_ITRF,TINYFRAMES_GCRF,TINYFRAMES_IAU2006A,matrix)==TINYFRAMES_OK);
    for(size_t i=0;i<3;++i){magnetic_ref[i]=0;for(size_t j=0;j<3;++j)magnetic_ref[i]+=matrix[i][j]*magnetic_ecef[j];}
    REQUIRE(tinysun_geocentric(tdb.first,tdb.fraction,&sun)==TINYSUN_OK);
    REQUIRE(tinysun_direction(&sun,rg,sun_ref)==TINYSUN_OK);
    /* Substitute real calibrated body-frame measurements in flight. */
    REQUIRE(tinyrotation_from_axis_angle(axis,.4,truth)==TINYROTATION_OK);
    REQUIRE(tinyrotation_apply(truth,magnetic_ref,magnetic_body)==TINYROTATION_OK);
    REQUIRE(tinyrotation_apply(truth,sun_ref,sun_body)==TINYROTATION_OK);
    REQUIRE(tinytriad_solve(magnetic_ref,sun_ref,magnetic_body,sun_body,attitude)==TINYTRIAD_OK);
    {double dot=0;for(size_t i=0;i<4;++i)dot+=truth[i]*attitude[i];REQUIRE(fabs(dot)>1-1e-12);}
    printf("TRIAD reference-to-body quaternion [w x y z]: %.9g %.9g %.9g %.9g\n",attitude[0],attitude[1],attitude[2],attitude[3]);
    return 0;
}
