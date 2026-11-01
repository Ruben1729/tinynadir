#include "tinycoverage/tinycoverage.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(x) do{if(!(x)){fprintf(stderr,"%s:%d %s\n",__FILE__,__LINE__,#x);exit(1);}}while(0)
#define PI 3.14159265358979323846
#include "data/reference.h"
static int access_margin(double t,void *ctx,double *out){(void)ctx;*out=cos(.001*t)-.2;return 0;}
static int broken(double t,void *ctx,double *out){(void)t;(void)ctx;(void)out;return -1;}
int main(void) {
 size_t i,n=999,calls=999;tinycoverage_arc arc;tinycoverage_point point;tinycoverage_footprint foot;tinycoverage_observation obs;tinycoverage_repeat repeat;
 tinycoverage_workspace scratch;tinycoverage_interval intervals[64],saved[64];
 for(i=0;i<sizeof arc_cases/sizeof arc_cases[0];i++) {
  CHECK(tinycoverage_great_circle(arc_cases[i][0],arc_cases[i][1],arc_cases[i][2],arc_cases[i][3],arc_cases[i][4],&arc)==TINYCOVERAGE_OK);
  CHECK(fabs(arc.surface_distance_m-arc_cases[i][5])<1e-7);CHECK(fabs(remainder(arc.forward_azimuth_rad-arc_cases[i][6],2*PI))<1e-12);
 }
 for(i=0;i<sizeof direct_cases/sizeof direct_cases[0];i++) {
  CHECK(tinycoverage_destination(direct_cases[i][0],direct_cases[i][1],direct_cases[i][2],direct_cases[i][3],direct_cases[i][4],&point)==TINYCOVERAGE_OK);
  CHECK(fabs(point.latitude_rad-direct_cases[i][5])<1e-12);CHECK(fabs(remainder(point.longitude_rad-direct_cases[i][6],2*PI))<1e-12);
 }
 for(i=0;i<sizeof footprint_cases/sizeof footprint_cases[0];i++) {
  CHECK(tinycoverage_nadir_footprint(footprint_cases[i][0],footprint_cases[i][1],footprint_cases[i][2],&foot)==TINYCOVERAGE_OK);
  CHECK(fabs(foot.central_angle_rad-footprint_cases[i][3])<1e-12);CHECK(fabs(foot.surface_radius_m-footprint_cases[i][4])<1e-5);CHECK(fabs(foot.area_m2-footprint_cases[i][5])<.1+1e-12*foot.area_m2);CHECK(foot.horizon_limited==(int)footprint_cases[i][6]);
 }
 for(i=0;i<sizeof observation_cases/sizeof observation_cases[0];i++) {
  CHECK(tinycoverage_observe(observation_cases[i],observation_cases[i][3],observation_cases[i][4],observation_cases[i]+5,&obs)==TINYCOVERAGE_OK);
  CHECK(fabs(obs.range_m-observation_cases[i][8])<1e-7);CHECK(fabs(remainder(obs.azimuth_rad-observation_cases[i][9],2*PI))<1e-12);CHECK(fabs(obs.elevation_rad-observation_cases[i][10])<1e-12);
 }
 CHECK(tinycoverage_repeat_orbit(3.986004418e14,86164.0905,15,1,&repeat)==TINYCOVERAGE_OK);CHECK(fabs(repeat.nodal_period_s*15-86164.0905)<1e-9);
 CHECK(tinycoverage_scan(access_margin,NULL,-1000,18000,300,1e-7,10000,&scratch,intervals,64,&n,&calls)==TINYCOVERAGE_OK);CHECK(n==sizeof access_cases/sizeof access_cases[0]);CHECK(calls<=10000);
 for(i=0;i<n;i++){CHECK(fabs(intervals[i].entry_s-access_cases[i][0])<1e-7);CHECK(fabs(intervals[i].exit_s-access_cases[i][1])<1e-7);}
 memcpy(saved,intervals,sizeof saved);n=999;calls=999;
 CHECK(tinycoverage_scan(access_margin,NULL,-1000,18000,300,1e-7,2,&scratch,intervals,64,&n,&calls)==TINYCOVERAGE_BUDGET);CHECK(n==999&&calls==999&&!memcmp(intervals,saved,sizeof saved));
 CHECK(tinycoverage_scan(access_margin,NULL,-1000,18000,300,1e-7,10000,&scratch,intervals,1,&n,&calls)==TINYCOVERAGE_CAPACITY);CHECK(n==999&&!memcmp(intervals,saved,sizeof saved));
 CHECK(tinycoverage_scan(broken,NULL,0,100,10,.01,100,&scratch,intervals,64,&n,&calls)==TINYCOVERAGE_CALLBACK);
 CHECK(tinycoverage_great_circle(6371e3,0,0,0,PI,&arc)==TINYCOVERAGE_GEOMETRY);
 CHECK(tinycoverage_observe(observation_cases[0],0,0,observation_cases[0],&obs)==TINYCOVERAGE_GEOMETRY);
 CHECK(tinycoverage_nadir_footprint(6371e3,-1,.1,&foot)==TINYCOVERAGE_RANGE);
 puts("tinycoverage: 264 independent geometry cases, SciPy access windows, and bounded-failure checks pass");return 0;
}
