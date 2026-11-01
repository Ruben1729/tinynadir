/* SPDX-License-Identifier: MIT */
#include "tinygeodesy/tinygeodesy.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "data/reference.h"
#define CHECK(x) do { if (!(x)) { fprintf(stderr,"line %d: %s\n",__LINE__,#x);return 1;} } while(0)
int main(void)
{
    size_t i,j;double max_pos=0,max_lat=0,max_height=0;
    for(i=0;i<sizeof reference/sizeof reference[0];++i) {
        tinygeodesy_ellipsoid e={reference[i][0],reference[i][1]};
        tinygeodesy_coordinate g={reference[i][2],reference[i][3],reference[i][4]},back;
        double r[3],C[9];
        CHECK(tinygeodesy_to_ecef(&e,&g,r)==TINYGEODESY_OK);
        CHECK(tinygeodesy_from_ecef(&e,r,&back)==TINYGEODESY_OK);
        for(j=0;j<3;++j){double error=fabs(r[j]-reference[i][5+j]);if(error>max_pos)max_pos=error;CHECK(error<2e-7);}
        if(fabs(back.latitude_rad-reference[i][8])>max_lat)max_lat=fabs(back.latitude_rad-reference[i][8]);
        if(fabs(back.height_m-reference[i][10])>max_height)max_height=fabs(back.height_m-reference[i][10]);
        CHECK(fabs(back.latitude_rad-g.latitude_rad)<2e-14);
        CHECK(fabs(back.longitude_rad-g.longitude_rad)<2e-14);
        CHECK(fabs(back.height_m-g.height_m)<3e-8);
        CHECK(fabs(back.latitude_rad-reference[i][8])<2e-14);CHECK(fabs(back.height_m-reference[i][10])<3e-8);
        CHECK(tinygeodesy_local_frame(g.latitude_rad,g.longitude_rad,TINYGEODESY_ENU,C)==TINYGEODESY_OK);
        for(j=0;j<9;++j)CHECK(fabs(C[j]-reference[i][11+j])<3e-15);
        CHECK(tinygeodesy_local_frame(g.latitude_rad,g.longitude_rad,TINYGEODESY_NED,C)==TINYGEODESY_OK);
        for(j=0;j<3;++j){CHECK(fabs(C[j]-reference[i][14+j])<3e-15);CHECK(fabs(C[3+j]-reference[i][11+j])<3e-15);CHECK(fabs(C[6+j]+reference[i][17+j])<3e-15);}
        CHECK(tinygeodesy_local_frame(g.latitude_rad,g.longitude_rad,TINYGEODESY_SEZ,C)==TINYGEODESY_OK);
        for(j=0;j<3;++j){CHECK(fabs(C[j]+reference[i][14+j])<3e-15);CHECK(fabs(C[3+j]-reference[i][11+j])<3e-15);CHECK(fabs(C[6+j]-reference[i][17+j])<3e-15);}
    }
    {
        tinygeodesy_ellipsoid e={TINYGEODESY_WGS84_A_M,TINYGEODESY_WGS84_F};
        tinygeodesy_coordinate g={1,2,3},saved=g;double r[3]={0,0,0},gc;
        CHECK(tinygeodesy_from_ecef(&e,r,&g)==TINYGEODESY_INVALID);CHECK(!memcmp(&g,&saved,sizeof g));
        r[2]=e.semi_major_axis_m*(1-e.flattening)+100;
        CHECK(tinygeodesy_from_ecef(&e,r,&g)==TINYGEODESY_OK);CHECK(g.longitude_rad==0);CHECK(fabs(g.height_m-100)<1e-8);
        CHECK(tinygeodesy_geocentric_latitude(.7,e.flattening,&gc)==TINYGEODESY_OK);
        CHECK(fabs(tan(gc)-(1-e.flattening)*(1-e.flattening)*tan(.7))<1e-15);
        CHECK(tinygeodesy_local_frame(NAN,0,TINYGEODESY_ENU,r)==TINYGEODESY_INVALID);
    }
    printf("Geodesy: 36 pymap3d cases; max position %.9g m, inverse comparison %.9g rad / %.9g m\n",max_pos,max_lat,max_height);
    return 0;
}
