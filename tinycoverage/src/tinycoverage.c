/* SPDX-License-Identifier: MIT; independent geometry/event-search implementation. */
#include "tinycoverage/tinycoverage.h"
#include <math.h>
#define PI 3.141592653589793238462643383279502884
static int angle(double x){return isfinite(x)&&fabs(x)<=PI/2;}
static double clamp(double x){return fmax(-1,fmin(1,x));}
static double azimuth(double x){x=fmod(x,2*PI);return x<0?x+2*PI:x;}
tinycoverage_status tinycoverage_great_circle(double radius,double lat1,double lon1,double lat2,double lon2,tinycoverage_arc *out) {
    double dlon,e,n,cross,c,central;tinycoverage_arc r;
    if(!out||!isfinite(radius)||!angle(lat1)||!angle(lat2)||!isfinite(lon1)||!isfinite(lon2))return TINYCOVERAGE_INVALID;
    if(radius<=0)return TINYCOVERAGE_RANGE;
    dlon=remainder(remainder(lon2,2*PI)-remainder(lon1,2*PI),2*PI);
    e=cos(lat2)*sin(dlon);n=cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(dlon);
    cross=hypot(e,n);c=sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(dlon);central=atan2(cross,c);
    if(cross<1e-14&&c<0)return TINYCOVERAGE_GEOMETRY;
    r.surface_distance_m=radius*central;r.forward_azimuth_rad=cross<1e-14?0:azimuth(atan2(e,n));
    if(!isfinite(r.surface_distance_m))return TINYCOVERAGE_RANGE;
    *out=r;return TINYCOVERAGE_OK;
}
tinycoverage_status tinycoverage_destination(double radius,double lat,double lon,double az,double range,tinycoverage_point *out) {
    double a,ca,sa,sl,cl,ta;tinycoverage_point r;
    if(!out||!isfinite(radius)||!angle(lat)||!isfinite(lon)||!isfinite(az)||!isfinite(range))return TINYCOVERAGE_INVALID;
    if(radius<=0||range<0||range/radius>PI)return TINYCOVERAGE_RANGE;
    a=range/radius;ca=cos(a);sa=sin(a);sl=sin(lat);cl=cos(lat);az=remainder(az,2*PI);
    /* Rodrigues rotation of the initial radial vector along the local tangent
       remains well-defined at a starting pole, with longitude defining north. */
    {double co=cos(remainder(lon,2*PI)),so=sin(remainder(lon,2*PI)),cz=cos(az),sz=sin(az);
     double x=ca*cl*co+sa*(-cz*sl*co-sz*so),y=ca*cl*so+sa*(-cz*sl*so+sz*co),z=ca*sl+sa*cz*cl;
     ta=hypot(x,y);r.latitude_rad=atan2(z,ta);r.longitude_rad=remainder(atan2(y,x),2*PI);}
    *out=r;return TINYCOVERAGE_OK;
}
tinycoverage_status tinycoverage_nadir_footprint(double radius,double altitude,double half,tinycoverage_footprint *out) {
    double satellite,horizon,limit,beta;tinycoverage_footprint r;
    if(!out||!isfinite(radius)||!isfinite(altitude)||!isfinite(half))return TINYCOVERAGE_INVALID;
    if(radius<=0||altitude<0||half<0||half>PI/2)return TINYCOVERAGE_RANGE;
    satellite=radius+altitude;if(!isfinite(satellite))return TINYCOVERAGE_RANGE;
    horizon=acos(clamp(radius/satellite));limit=asin(clamp(radius/satellite));
    r.horizon_limited=half>=limit;
    beta=r.horizon_limited?horizon:asin(clamp((satellite/radius)*sin(half)))-half;
    r.central_angle_rad=fmax(0,beta);r.surface_radius_m=radius*r.central_angle_rad;r.area_m2=4*PI*radius*radius*pow(sin(r.central_angle_rad/2),2);
    if(!isfinite(r.area_m2)||!isfinite(r.surface_radius_m))return TINYCOVERAGE_RANGE;
    *out=r;return TINYCOVERAGE_OK;
}
tinycoverage_status tinycoverage_observe(const double site[3],double lat,double lon,const double object[3],tinycoverage_observation *out) {
    double v[3],range,co,so,cl,sl,e,n,u;tinycoverage_observation r;int i;
    if(!site||!object||!out||!angle(lat)||!isfinite(lon))return TINYCOVERAGE_INVALID;
    for(i=0;i<3;i++){if(!isfinite(site[i])||!isfinite(object[i]))return TINYCOVERAGE_INVALID;v[i]=object[i]-site[i];}
    range=hypot(hypot(v[0],v[1]),v[2]);if(!isfinite(range))return TINYCOVERAGE_RANGE;if(range==0)return TINYCOVERAGE_GEOMETRY;
    lon=remainder(lon,2*PI);co=cos(lon);so=sin(lon);cl=cos(lat);sl=sin(lat);
    e=-so*v[0]+co*v[1];n=-sl*co*v[0]-sl*so*v[1]+cl*v[2];u=cl*co*v[0]+cl*so*v[1]+sl*v[2];
    r.range_m=range;r.azimuth_rad=hypot(e,n)<1e-12*range?0:azimuth(atan2(e,n));r.elevation_rad=atan2(u,hypot(e,n));
    *out=r;return TINYCOVERAGE_OK;
}
tinycoverage_status tinycoverage_repeat_orbit(double mu,double day,unsigned rev,unsigned turns,tinycoverage_repeat *out) {
    tinycoverage_repeat r;double n;
    if(!out||!isfinite(mu)||!isfinite(day))return TINYCOVERAGE_INVALID;
    if(mu<=0||day<=0||rev<1||turns<1||rev>100000||turns>100000)return TINYCOVERAGE_RANGE;
    r.nodal_period_s=day*((double)turns/rev);n=2*PI/r.nodal_period_s;r.semimajor_axis_m=cbrt(mu)/pow(n,2./3.);
    if(!isfinite(r.nodal_period_s)||r.nodal_period_s<=0||!isfinite(r.semimajor_axis_m)||r.semimajor_axis_m<=0)return TINYCOVERAGE_RANGE;
    *out=r;return TINYCOVERAGE_OK;
}
static tinycoverage_status evaluate(tinycoverage_margin_fn fn,void *context,double t,size_t maximum,size_t *calls,double *value) {
    if(*calls>=maximum)return TINYCOVERAGE_BUDGET;
    ++*calls;if(fn(t,context,value)||!isfinite(*value))return TINYCOVERAGE_CALLBACK;
    return TINYCOVERAGE_OK;
}
static tinycoverage_status crossing(tinycoverage_margin_fn fn,void *ctx,double a,double b,double fa,double fb,double tolerance,size_t budget,size_t *calls,double *event) {
    unsigned k;tinycoverage_status status;int side=fa>=0;
    if(fa==0){*event=a;return TINYCOVERAGE_OK;}if(fb==0){*event=b;return TINYCOVERAGE_OK;}
    for(k=0;k<64;k++) {
        double m,f;if(b-a<=tolerance){*event=a+(b-a)/2;return TINYCOVERAGE_OK;}
        m=a+(b-a)/2;if(m==a||m==b)return TINYCOVERAGE_RESOLUTION;
        status=evaluate(fn,ctx,m,budget,calls,&f);if(status!=TINYCOVERAGE_OK)return status;
        if(f==0){*event=m;return TINYCOVERAGE_OK;}if((f>=0)==side)a=m;else b=m;
    }
    return TINYCOVERAGE_RESOLUTION;
}
tinycoverage_status tinycoverage_scan(tinycoverage_margin_fn fn,void *ctx,double start,double end,double sample,double tolerance,size_t budget,tinycoverage_workspace *scratch,tinycoverage_interval *out,size_t capacity,size_t *count,size_t *calls) {
    size_t used=0,invocations=0,i;double t,f,opened=0;tinycoverage_status status;int inside;
    if(!fn||!scratch||!out||!count||!calls||!isfinite(start)||!isfinite(end)||!isfinite(sample)||!isfinite(tolerance))return TINYCOVERAGE_INVALID;
    if(end<start||sample<=0||tolerance<=0||capacity<1||capacity>TINYCOVERAGE_MAX_INTERVALS||budget<1||budget>1000000)return TINYCOVERAGE_RANGE;
    t=start;status=evaluate(fn,ctx,t,budget,&invocations,&f);if(status!=TINYCOVERAGE_OK)return status;
    inside=f>=0;if(inside)opened=t;
    while(t<end) {
        double next=fmin(end,t+sample),value,event;int entering;
        if(!isfinite(next)||next<=t)return TINYCOVERAGE_RESOLUTION;
        status=evaluate(fn,ctx,next,budget,&invocations,&value);if(status!=TINYCOVERAGE_OK)return status;
        entering=value>=0;
        if(entering!=inside) {
            status=crossing(fn,ctx,t,next,f,value,tolerance,budget,&invocations,&event);if(status!=TINYCOVERAGE_OK)return status;
            if(entering)opened=event;
            else {if(used>=capacity)return TINYCOVERAGE_CAPACITY;scratch->intervals[used].entry_s=opened;scratch->intervals[used++].exit_s=event;}
        }
        inside=entering;t=next;f=value;
    }
    if(inside){if(used>=capacity)return TINYCOVERAGE_CAPACITY;scratch->intervals[used].entry_s=opened;scratch->intervals[used++].exit_s=end;}
    for(i=0;i<used;i++)out[i]=scratch->intervals[i];
    *count=used;*calls=invocations;return TINYCOVERAGE_OK;
}
