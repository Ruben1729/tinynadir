/* SPDX-License-Identifier: MIT */
#include "tinyobservations/tinyobservations.h"
#include <math.h>
#define TO_PI 3.1415926535897932384626433832795
static int valid_state(const tinyobservations_state *s)
{unsigned i;if(!s)return 0;for(i=0;i<3;++i)if(!isfinite(s->position_m[i])||!isfinite(s->velocity_m_s[i]))return 0;return 1;}
tinyobservations_status tinyobservations_to_radec(const tinyobservations_state *s,tinyobservations_radec *out)
{
    tinyobservations_radec r;double x,y,z,p,radius,xd,yd,zd,pdot;
    if(!out || !valid_state(s))return TINYOBSERVATIONS_INVALID;
    x=s->position_m[0];y=s->position_m[1];z=s->position_m[2];xd=s->velocity_m_s[0];yd=s->velocity_m_s[1];zd=s->velocity_m_s[2];
    p=hypot(x,y);radius=hypot(p,z);
    if(!(radius>0) || !isfinite(radius))return TINYOBSERVATIONS_INVALID;
    if(p/radius<1e-14)return TINYOBSERVATIONS_SINGULAR;
    r.range_m=radius;r.right_ascension_rad=atan2(y,x);if(r.right_ascension_rad<0)r.right_ascension_rad+=2*TO_PI;
    r.declination_rad=atan2(z,p);
    r.range_rate_m_s=(x/radius)*xd+(y/radius)*yd+(z/radius)*zd;
    pdot=(x/p)*xd+(y/p)*yd;
    r.right_ascension_rate_rad_s=((x/p)*yd-(y/p)*xd)/p;
    r.declination_rate_rad_s=((p/radius)*zd-(z/radius)*pdot)/radius;
    if(!isfinite(r.range_rate_m_s) || !isfinite(r.right_ascension_rate_rad_s) || !isfinite(r.declination_rate_rad_s))return TINYOBSERVATIONS_INVALID;
    *out=r;return TINYOBSERVATIONS_OK;
}
tinyobservations_status tinyobservations_from_radec(const tinyobservations_radec *r,tinyobservations_state *out)
{
    tinyobservations_state s;double a,d,ca,sa,cd,sd;unsigned i;
    if(!r || !out || !isfinite(r->range_m) || r->range_m<=0 || !isfinite(r->right_ascension_rad) ||
       !isfinite(r->declination_rad) || fabs(r->declination_rad)>TO_PI/2 || !isfinite(r->range_rate_m_s) ||
       !isfinite(r->right_ascension_rate_rad_s) || !isfinite(r->declination_rate_rad_s))return TINYOBSERVATIONS_INVALID;
    a=r->right_ascension_rad;d=r->declination_rad;ca=cos(a);sa=sin(a);cd=cos(d);sd=sin(d);
    s.position_m[0]=r->range_m*cd*ca;s.position_m[1]=r->range_m*cd*sa;s.position_m[2]=r->range_m*sd;
    s.velocity_m_s[0]=r->range_rate_m_s*cd*ca-r->range_m*(sd*ca*r->declination_rate_rad_s+cd*sa*r->right_ascension_rate_rad_s);
    s.velocity_m_s[1]=r->range_rate_m_s*cd*sa+r->range_m*(-sd*sa*r->declination_rate_rad_s+cd*ca*r->right_ascension_rate_rad_s);
    s.velocity_m_s[2]=r->range_rate_m_s*sd+r->range_m*cd*r->declination_rate_rad_s;
    for(i=0;i<3;++i)if(!isfinite(s.position_m[i])||!isfinite(s.velocity_m_s[i]))return TINYOBSERVATIONS_INVALID;
    *out=s;return TINYOBSERVATIONS_OK;
}
tinyobservations_status tinyobservations_to_azel(const tinygeodesy_ellipsoid *e,const tinygeodesy_coordinate *station,const tinyobservations_state *s,tinyobservations_azel *out)
{
    double origin[3],C[9],enu[3]={0},enud[3]={0};tinyobservations_state local;tinyobservations_radec a;tinyobservations_azel result;unsigned i,j;tinyobservations_status status;
    if(!out || !valid_state(s) || tinygeodesy_to_ecef(e,station,origin)!=TINYGEODESY_OK || tinygeodesy_local_frame(station->latitude_rad,station->longitude_rad,TINYGEODESY_ENU,C)!=TINYGEODESY_OK)return TINYOBSERVATIONS_INVALID;
    for(i=0;i<3;++i)for(j=0;j<3;++j){enu[i]+=C[3*i+j]*(s->position_m[j]-origin[j]);enud[i]+=C[3*i+j]*s->velocity_m_s[j];}
    local.position_m[0]=enu[1];local.position_m[1]=enu[0];local.position_m[2]=enu[2];local.velocity_m_s[0]=enud[1];local.velocity_m_s[1]=enud[0];local.velocity_m_s[2]=enud[2];
    status=tinyobservations_to_radec(&local,&a);if(status!=TINYOBSERVATIONS_OK)return status;
    result.range_m=a.range_m;result.azimuth_rad=a.right_ascension_rad;result.elevation_rad=a.declination_rad;result.range_rate_m_s=a.range_rate_m_s;result.azimuth_rate_rad_s=a.right_ascension_rate_rad_s;result.elevation_rate_rad_s=a.declination_rate_rad_s;
    *out=result;return TINYOBSERVATIONS_OK;
}
tinyobservations_status tinyobservations_from_azel(const tinygeodesy_ellipsoid *e,const tinygeodesy_coordinate *station,const tinyobservations_azel *angles,tinyobservations_state *out)
{
    tinyobservations_radec a;tinyobservations_state local,result;double origin[3],C[9],enu[3],enud[3];unsigned i,j;tinyobservations_status status;
    if(!angles || !out || tinygeodesy_to_ecef(e,station,origin)!=TINYGEODESY_OK || tinygeodesy_local_frame(station->latitude_rad,station->longitude_rad,TINYGEODESY_ENU,C)!=TINYGEODESY_OK)return TINYOBSERVATIONS_INVALID;
    a.range_m=angles->range_m;a.right_ascension_rad=angles->azimuth_rad;a.declination_rad=angles->elevation_rad;a.range_rate_m_s=angles->range_rate_m_s;a.right_ascension_rate_rad_s=angles->azimuth_rate_rad_s;a.declination_rate_rad_s=angles->elevation_rate_rad_s;
    status=tinyobservations_from_radec(&a,&local);if(status!=TINYOBSERVATIONS_OK)return status;
    enu[0]=local.position_m[1];enu[1]=local.position_m[0];enu[2]=local.position_m[2];enud[0]=local.velocity_m_s[1];enud[1]=local.velocity_m_s[0];enud[2]=local.velocity_m_s[2];
    for(i=0;i<3;++i){result.position_m[i]=origin[i];result.velocity_m_s[i]=0;for(j=0;j<3;++j){result.position_m[i]+=C[3*j+i]*enu[j];result.velocity_m_s[i]+=C[3*j+i]*enud[j];}}
    if(!valid_state(&result))return TINYOBSERVATIONS_INVALID;
    *out=result;return TINYOBSERVATIONS_OK;
}
tinyobservations_status tinyobservations_relative(const tinyobservations_state *s,const tinyobservations_state *observer,tinyobservations_state *out)
{
    tinyobservations_state r;unsigned i;if(!out||!valid_state(s)||!valid_state(observer))return TINYOBSERVATIONS_INVALID;
    for(i=0;i<3;++i){r.position_m[i]=s->position_m[i]-observer->position_m[i];r.velocity_m_s[i]=s->velocity_m_s[i]-observer->velocity_m_s[i];}
    if(!valid_state(&r))return TINYOBSERVATIONS_INVALID;*out=r;return TINYOBSERVATIONS_OK;
}
tinyobservations_status tinyobservations_ecliptic(const tinyobservations_state *s,double eps,tinyobservations_state *out)
{
    tinyobservations_state r;double c,d;if(!out||!valid_state(s)||!isfinite(eps))return TINYOBSERVATIONS_INVALID;c=cos(eps);d=sin(eps);
    r.position_m[0]=s->position_m[0];r.position_m[1]=c*s->position_m[1]+d*s->position_m[2];r.position_m[2]=-d*s->position_m[1]+c*s->position_m[2];
    r.velocity_m_s[0]=s->velocity_m_s[0];r.velocity_m_s[1]=c*s->velocity_m_s[1]+d*s->velocity_m_s[2];r.velocity_m_s[2]=-d*s->velocity_m_s[1]+c*s->velocity_m_s[2];
    if(!valid_state(&r))return TINYOBSERVATIONS_INVALID;*out=r;return TINYOBSERVATIONS_OK;
}
