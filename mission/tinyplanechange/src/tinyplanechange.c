/* SPDX-License-Identifier: MIT */
#include "tinyplanechange/tinyplanechange.h"
#include <math.h>
#define TP_PI 3.1415926535897932384626433832795
static double norm(const double x[3]){return hypot(hypot(x[0],x[1]),x[2]);}
static int valid(double v,double fpa){return isfinite(v)&&v>=0&&isfinite(fpa)&&fabs(fpa)<=TP_PI/2;}
static double wrap(double x){x=fmod(x,2*TP_PI);return x<0 ? x+2*TP_PI : x;}
tinyplanechange_status tinyplanechange_inclination(double v,double fpa,double di,double *out)
{
    double dv;
    if(!out || !valid(v,fpa) || !isfinite(di))return TINYPLANECHANGE_INVALID;
    dv=(2*cos(fpa)*fabs(sin(di/2)))*v;
    if(!isfinite(dv))return TINYPLANECHANGE_INVALID;
    *out=dv;return TINYPLANECHANGE_OK;
}
tinyplanechange_status tinyplanechange_combined(double v,double fpa,double i1,double i2,double O,tinyplanechange_result *out)
{
    double h1[3],h2[3],plus[3],minus[3],r[3],rn;unsigned j;tinyplanechange_result result;
    if(!out || !valid(v,fpa) || !isfinite(i1) || i1<0 || i1>TP_PI || !isfinite(i2) || i2<0 || i2>TP_PI || !isfinite(O))return TINYPLANECHANGE_INVALID;
    h1[0]=0;h1[1]=-sin(i1);h1[2]=cos(i1);h2[0]=sin(i2)*sin(O);h2[1]=-sin(i2)*cos(O);h2[2]=cos(i2);
    for(j=0;j<3;++j){plus[j]=h1[j]+h2[j];minus[j]=h1[j]-h2[j];}
    result.plane_angle_rad=2*atan2(norm(minus),norm(plus));
    result.delta_v_m_s=(2*cos(fpa)*sin(result.plane_angle_rad/2))*v;
    if(!isfinite(result.delta_v_m_s))return TINYPLANECHANGE_INVALID;
    r[0]=h1[1]*h2[2]-h1[2]*h2[1];r[1]=h1[2]*h2[0]-h1[0]*h2[2];r[2]=h1[0]*h2[1]-h1[1]*h2[0];rn=norm(r);
    result.unique_intersection=rn>1e-12;
    if(result.unique_intersection) {
        result.initial_argument_latitude_rad=wrap(atan2(r[1]*cos(i1)+r[2]*sin(i1),r[0]));
        result.final_argument_latitude_rad=wrap(atan2(-r[0]*cos(i2)*sin(O)+r[1]*cos(i2)*cos(O)+r[2]*sin(i2),r[0]*cos(O)+r[1]*sin(O)));
    } else {result.initial_argument_latitude_rad=0;result.final_argument_latitude_rad=0;}
    *out=result;return TINYPLANECHANGE_OK;
}
tinyplanechange_status tinyplanechange_node(double v,double fpa,double i,double O,tinyplanechange_result *out)
{return tinyplanechange_combined(v,fpa,i,i,O,out);}
