/* SPDX-License-Identifier: MIT */
#include "tinyobservations/tinyobservations.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "data/reference.h"
#define CHECK(x) do { if (!(x)) {fprintf(stderr,"line %d: %s\n",__LINE__,#x);return 1;} } while(0)
int main(void)
{
    size_t i,j;double max_range=0,max_angle=0,max_rate=0;
    for(i=0;i<sizeof radec_reference/sizeof radec_reference[0];++i){const double *r=radec_reference[i];tinyobservations_state s,back,ecl;tinyobservations_radec a;for(j=0;j<3;++j){s.position_m[j]=r[j];s.velocity_m_s[j]=r[3+j];}CHECK(tinyobservations_to_radec(&s,&a)==TINYOBSERVATIONS_OK);CHECK(fabs(a.range_m-r[6])<1e-7);CHECK(fabs(remainder(a.right_ascension_rad-r[7],6.2831853071795864769))<2e-14);CHECK(fabs(a.declination_rad-r[8])<2e-14);CHECK(fabs(a.range_rate_m_s-r[9])<1e-8);CHECK(fabs(a.right_ascension_rate_rad_s-r[10])<1e-14);CHECK(fabs(a.declination_rate_rad_s-r[11])<1e-14);CHECK(tinyobservations_from_radec(&a,&back)==TINYOBSERVATIONS_OK);for(j=0;j<3;++j){CHECK(fabs(back.position_m[j]-s.position_m[j])<1e-7);CHECK(fabs(back.velocity_m_s[j]-s.velocity_m_s[j])<1e-9);}CHECK(tinyobservations_ecliptic(&s,.409,&ecl)==TINYOBSERVATIONS_OK);CHECK(tinyobservations_ecliptic(&ecl,-.409,&back)==TINYOBSERVATIONS_OK);for(j=0;j<3;++j){CHECK(fabs(back.position_m[j]-s.position_m[j])<1e-7);CHECK(fabs(back.velocity_m_s[j]-s.velocity_m_s[j])<1e-9);}}
    for(i=0;i<sizeof azel_reference/sizeof azel_reference[0];++i){const double *r=azel_reference[i];tinygeodesy_ellipsoid ell={TINYGEODESY_WGS84_A_M,TINYGEODESY_WGS84_F};tinygeodesy_coordinate station={r[0],r[1],r[2]};tinyobservations_state s,back;tinyobservations_azel a;double er,ea,ev;for(j=0;j<3;++j){s.position_m[j]=r[3+j];s.velocity_m_s[j]=r[6+j];}CHECK(tinyobservations_to_azel(&ell,&station,&s,&a)==TINYOBSERVATIONS_OK);er=fabs(a.range_m-r[9]);ea=fmax(fabs(remainder(a.azimuth_rad-r[10],6.2831853071795864769)),fabs(a.elevation_rad-r[11]));ev=fabs(a.range_rate_m_s-r[12]);if(er>max_range)max_range=er;if(ea>max_angle)max_angle=ea;if(ev>max_rate)max_rate=ev;CHECK(er<1e-7);CHECK(ea<2e-14);CHECK(ev<1e-8);CHECK(fabs(a.azimuth_rate_rad_s-r[13])<1e-14);CHECK(fabs(a.elevation_rate_rad_s-r[14])<1e-14);CHECK(tinyobservations_from_azel(&ell,&station,&a,&back)==TINYOBSERVATIONS_OK);for(j=0;j<3;++j){CHECK(fabs(back.position_m[j]-s.position_m[j])<1e-7);CHECK(fabs(back.velocity_m_s[j]-s.velocity_m_s[j])<1e-9);}}
    {tinyobservations_state s={{0,0,1},{0,0,0}};tinyobservations_radec a={1,2,3,4,5,6},saved=a;CHECK(tinyobservations_to_radec(&s,&a)==TINYOBSERVATIONS_SINGULAR);CHECK(!memcmp(&a,&saved,sizeof a));s.position_m[2]=0;CHECK(tinyobservations_to_radec(&s,&a)==TINYOBSERVATIONS_INVALID);CHECK(!memcmp(&a,&saved,sizeof a));}
    printf("Observations: 3 Skyfield RA/Dec +9 pymap3d az/el cases; max %.9g m / %.9g rad / %.9g m/s\n",max_range,max_angle,max_rate);return 0;
}
