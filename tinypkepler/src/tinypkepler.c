/* SPDX-License-Identifier: MIT */
#include "tinypkepler/tinypkepler.h"
#include "tinyanomaly/tinyanomaly.h"
#include <math.h>
#define TK_PI 3.1415926535897932384626433832795
static double wrap(double x){x=fmod(x,2*TK_PI);return x<0 ? x+2*TK_PI : x;}
tinypkepler_status tinypkepler_propagate(double mu,double re,double j2,const tinyelements_classical *initial,double dt,tinypkepler_result *out)
{
    tinypkepler_result result;tinyelements_state check;tinyelements_classical e;
    tinyanomaly_solution phase;double a,n,k,M,E,c,ratio;
    if(!out||!initial||!isfinite(mu)||mu<=0||!isfinite(re)||re<=0||!isfinite(j2)||!isfinite(dt)||
       initial->eccentricity>=1||tinyelements_from_classical(mu,initial,&check)!=TINYELEMENTS_OK)return TINYPKEPLER_INVALID;
    e=*initial;a=e.p_m/((1-e.eccentricity)*(1+e.eccentricity));n=sqrt(mu/a)/a;
    if(!isfinite(a)||!(a>0)||!isfinite(n)||!(n>0)||!isfinite(n*dt)||fabs(n*dt)>1e8)return TINYPKEPLER_INVALID;
    ratio=re/e.p_m;k=1.5*j2*n*ratio*ratio;c=cos(e.inclination_rad);
    result.raan_rate_rad_s=-k*c;result.arg_periapsis_rate_rad_s=.5*k*(5*c*c-1);result.mean_motion_rad_s=n;
    if(!isfinite(result.raan_rate_rad_s)||!isfinite(result.arg_periapsis_rate_rad_s))return TINYPKEPLER_INVALID;
    if(tinyanomaly_from_true(e.eccentricity,remainder(e.true_anomaly_rad,2*TK_PI),&M,&E)!=TINYANOMALY_OK)return TINYPKEPLER_INVALID;
    if(tinyanomaly_from_mean(e.eccentricity,M+n*dt,&phase)!=TINYANOMALY_OK)return TINYPKEPLER_NO_CONVERGENCE;
    e.raan_rad=wrap(e.raan_rad+result.raan_rate_rad_s*dt);
    e.arg_periapsis_rad=wrap(e.arg_periapsis_rad+result.arg_periapsis_rate_rad_s*dt);
    e.true_anomaly_rad=wrap(phase.true_anomaly_rad);
    if(tinyelements_from_classical(mu,&e,&result.state)!=TINYELEMENTS_OK)return TINYPKEPLER_INVALID;
    result.mean_elements=e;*out=result;return TINYPKEPLER_OK;
}
