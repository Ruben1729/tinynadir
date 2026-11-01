/* SPDX-License-Identifier: MIT */
#include "tinyhohmann/tinyhohmann.h"
#include "tinyanomaly/tinyanomaly.h"
#include <math.h>
#define TH_PI 3.1415926535897932384626433832795
static int valid(double mu,double r){return isfinite(mu)&&mu>0&&isfinite(r)&&r>0;}
static double half_period(double mu,double a){return TH_PI*a*sqrt(a/mu);}
static int finite_transfer(const tinyhohmann_transfer *t)
{return isfinite(t->burn_m_s[0])&&isfinite(t->burn_m_s[1])&&isfinite(t->burn_m_s[2])&&isfinite(t->total_delta_v_m_s)&&isfinite(t->time_of_flight_s)&&t->time_of_flight_s>=0;}
tinyhohmann_status tinyhohmann_circular(double mu,double r1,double r2,tinyhohmann_transfer *out)
{
    tinyhohmann_transfer t;double a,v1,v2;
    if(!out || !valid(mu,r1) || !valid(mu,r2))return TINYHOHMANN_INVALID;
    a=.5*r1+.5*r2;v1=sqrt(mu/r1);v2=sqrt(mu/r2);
    t.burn_m_s[0]=fabs(v1*(sqrt(r2/a)-1));t.burn_m_s[1]=fabs(v2*(1-sqrt(r1/a)));t.burn_m_s[2]=0;
    t.total_delta_v_m_s=t.burn_m_s[0]+t.burn_m_s[1];t.time_of_flight_s=half_period(mu,a);t.burn_count=2;
    if(!finite_transfer(&t))return TINYHOHMANN_INVALID;
    *out=t;return TINYHOHMANN_OK;
}
tinyhohmann_status tinyhohmann_bielliptic(double mu,double r1,double rb,double r2,tinyhohmann_transfer *out)
{
    tinyhohmann_transfer t;double a1,a2,vb;
    if(!out || !valid(mu,r1) || !valid(mu,r2) || !valid(mu,rb))return TINYHOHMANN_INVALID;
    if(rb<r1 || rb<r2)return TINYHOHMANN_GEOMETRY;
    a1=.5*r1+.5*rb;a2=.5*r2+.5*rb;vb=sqrt(mu/rb);
    t.burn_m_s[0]=fabs(sqrt(mu/r1)*(sqrt(rb/a1)-1));
    t.burn_m_s[1]=fabs(vb*(sqrt(r2/a2)-sqrt(r1/a1)));
    t.burn_m_s[2]=fabs(sqrt(mu/r2)*(1-sqrt(rb/a2)));
    t.total_delta_v_m_s=t.burn_m_s[0]+t.burn_m_s[1]+t.burn_m_s[2];
    t.time_of_flight_s=half_period(mu,a1)+half_period(mu,a2);t.burn_count=3;
    if(!finite_transfer(&t))return TINYHOHMANN_INVALID;
    *out=t;return TINYHOHMANN_OK;
}
tinyhohmann_status tinyhohmann_one_tangent(double mu,double r1,double r2,int apo,double nu,tinyhohmann_tangent *out)
{
    tinyhohmann_tangent t;double c0=apo ? -1.0 : 1.0,c=cos(nu),den,e,p,vscale,vt1,vt2,vr2,M,A,dt;
    if(!out || !valid(mu,r1) || !valid(mu,r2) || !isfinite(nu) ||
       (!apo && !(nu>0 && nu<=TH_PI)) || (apo && !(nu>TH_PI && nu<2*TH_PI)))return TINYHOHMANN_INVALID;
    den=(r2/r1)*c-c0;
    if(den==0)return TINYHOHMANN_GEOMETRY;
    e=(1-r2/r1)/den;p=r1*(1+e*c0);
    if(!isfinite(e) || e<0 || !(p>0) || !(1+e*c>0) || !isfinite(p))return TINYHOHMANN_GEOMETRY;
    if(tinyanomaly_from_true(e,nu>TH_PI ? nu-2*TH_PI : nu,&M,&A)!=TINYANOMALY_OK)return TINYHOHMANN_GEOMETRY;
    if(e<1) {
        double a=p/((1-e)*(1+e)),M0=apo ? TH_PI : 0;
        if(M<0)M+=2*TH_PI;
        dt=(M-M0)*a*sqrt(a/mu);
    } else if(e==1)dt=.5*p*sqrt(p/mu)*M;
    else {double aa=p/((e-1)*(e+1));dt=M*aa*sqrt(aa/mu);}
    if(!(dt>0) || !isfinite(dt))return TINYHOHMANN_GEOMETRY;
    vscale=sqrt(mu/p);vt1=vscale*(1+e*c0);vt2=vscale*(1+e*c);vr2=vscale*e*sin(nu);
    t.transfer.burn_m_s[0]=fabs(vt1-sqrt(mu/r1));
    t.transfer.burn_m_s[1]=hypot(vr2,sqrt(mu/r2)-vt2);t.transfer.burn_m_s[2]=0;
    t.transfer.total_delta_v_m_s=t.transfer.burn_m_s[0]+t.transfer.burn_m_s[1];t.transfer.burn_count=2;t.transfer.time_of_flight_s=dt;
    t.eccentricity=e;t.p_m=p;t.arrival_flight_path_rad=atan2(vr2,vt2);
    if(!finite_transfer(&t.transfer))return TINYHOHMANN_INVALID;
    *out=t;return TINYHOHMANN_OK;
}
