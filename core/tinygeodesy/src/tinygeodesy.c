/* SPDX-License-Identifier: MIT */
#include "tinygeodesy/tinygeodesy.h"
#include <math.h>
#include <float.h>
#define TG_PI 3.1415926535897932384626433832795
static int valid(const tinygeodesy_ellipsoid *e)
{ return e && isfinite(e->semi_major_axis_m) && e->semi_major_axis_m>0.0 &&
    isfinite(e->flattening) && e->flattening>=0.0 && e->flattening<=0.1; }
tinygeodesy_status tinygeodesy_to_ecef(const tinygeodesy_ellipsoid *e,
                                     const tinygeodesy_coordinate *g,double out[3])
{
    double s,c,e2,N,v[3];unsigned i;
    if (!valid(e) || !g || !out || !isfinite(g->latitude_rad) || fabs(g->latitude_rad)>TG_PI/2 ||
        !isfinite(g->longitude_rad) || !isfinite(g->height_m) || g->height_m< -e->semi_major_axis_m/2)
        return TINYGEODESY_INVALID;
    s=sin(g->latitude_rad);c=cos(g->latitude_rad);e2=e->flattening*(2-e->flattening);
    N=e->semi_major_axis_m/sqrt(1-e2*s*s);
    v[0]=(N+g->height_m)*c*cos(g->longitude_rad);
    v[1]=(N+g->height_m)*c*sin(g->longitude_rad);
    v[2]=(N*(1-e2)+g->height_m)*s;
    for(i=0;i<3;++i) if (!isfinite(v[i])) return TINYGEODESY_INVALID;
    for(i=0;i<3;++i) out[i]=v[i];
    return TINYGEODESY_OK;
}
tinygeodesy_status tinygeodesy_from_ecef(const tinygeodesy_ellipsoid *e,
                                       const double r[3],tinygeodesy_coordinate *out)
{
    double p,norm,e2,lat;unsigned iteration;
    tinygeodesy_coordinate result;
    if (!valid(e) || !r || !out || !isfinite(r[0]) || !isfinite(r[1]) || !isfinite(r[2]))
        return TINYGEODESY_INVALID;
    p=hypot(r[0],r[1]);norm=hypot(p,r[2]);
    if (!isfinite(norm) || norm<e->semi_major_axis_m/4) return TINYGEODESY_INVALID;
    e2=e->flattening*(2-e->flattening);
    lat=atan2(r[2],p*(1-e2));
    for(iteration=0;iteration<100;++iteration) {
        double s=sin(lat),N=e->semi_major_axis_m/sqrt(1-e2*s*s);
        double next=atan2(r[2]+e2*N*s,p);
        if(fabs(next-lat)<=2*DBL_EPSILON) {lat=next;break;}
        lat=next;
    }
    if(iteration==100) return TINYGEODESY_NO_CONVERGENCE;
    result.latitude_rad=lat;result.longitude_rad=p==0.0 ? 0.0 : atan2(r[1],r[0]);
    result.height_m=p*cos(lat)+r[2]*sin(lat)-e->semi_major_axis_m*sqrt(1-e2*sin(lat)*sin(lat));
    if(!isfinite(result.height_m)) return TINYGEODESY_INVALID;
    *out=result;return TINYGEODESY_OK;
}
tinygeodesy_status tinygeodesy_local_frame(double lat,double lon,tinygeodesy_frame frame,double out[9])
{
    double E[3],N[3],U[3],result[9],s,c,sl,cl;unsigned i;
    if(!out || !isfinite(lat) || fabs(lat)>TG_PI/2 || !isfinite(lon) ||
       (frame!=TINYGEODESY_ENU && frame!=TINYGEODESY_NED && frame!=TINYGEODESY_SEZ)) return TINYGEODESY_INVALID;
    s=sin(lat);c=cos(lat);sl=sin(lon);cl=cos(lon);
    E[0]=-sl;E[1]=cl;E[2]=0;N[0]=-s*cl;N[1]=-s*sl;N[2]=c;U[0]=c*cl;U[1]=c*sl;U[2]=s;
    for(i=0;i<3;++i) {
        if(frame==TINYGEODESY_ENU) {result[i]=E[i];result[3+i]=N[i];result[6+i]=U[i];}
        else if(frame==TINYGEODESY_NED) {result[i]=N[i];result[3+i]=E[i];result[6+i]=-U[i];}
        else {result[i]=-N[i];result[3+i]=E[i];result[6+i]=U[i];}
    }
    for(i=0;i<9;++i) out[i]=result[i];
    return TINYGEODESY_OK;
}
tinygeodesy_status tinygeodesy_geocentric_latitude(double lat,double f,double *out)
{
    if(!out || !isfinite(lat) || fabs(lat)>TG_PI/2 || !isfinite(f) || f<0 || f>0.1) return TINYGEODESY_INVALID;
    *out=atan2((1-f)*(1-f)*sin(lat),cos(lat));
    return TINYGEODESY_OK;
}
