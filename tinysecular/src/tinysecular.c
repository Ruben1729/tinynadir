/* SPDX-License-Identifier: MIT; independently written first-order J2 equations. */
#include "tinysecular/tinysecular.h"
#include <math.h>
#define PI 3.141592653589793238462643383279502884
static double wrap(double x){double r=fmod(x,2*PI);return r<0?r+2*PI:r;}
tinysecular_status tinysecular_j2_rates(const tinysecular_model *m,const tinysecular_elements *e,tinysecular_profile profile,tinysecular_rates *out) {
 double a,ecc,c,n,p,factor; tinysecular_rates r;
 if(!m||!e||!out)return TINYSECULAR_INVALID;
 if(!isfinite(m->mu_m3_s2)||!isfinite(m->equatorial_radius_m)||!isfinite(m->j2)||!isfinite(e->semimajor_axis_m)||!isfinite(e->eccentricity)||!isfinite(e->inclination_rad)||!isfinite(e->raan_rad)||!isfinite(e->arg_periapsis_rad)||!isfinite(e->mean_anomaly_rad))return TINYSECULAR_INVALID;
 a=e->semimajor_axis_m;ecc=e->eccentricity;
 if(m->mu_m3_s2<=0||m->equatorial_radius_m<=0||fabs(m->j2)>.01||a<=0||ecc<0||ecc>=1||e->inclination_rad<0||e->inclination_rad>PI||a*(1-ecc)<=m->equatorial_radius_m||fabs(e->raan_rad)>1e12||fabs(e->arg_periapsis_rad)>1e12||fabs(e->mean_anomaly_rad)>1e12||(profile!=TINYSECULAR_AVERAGED_J2&&profile!=TINYSECULAR_VALLADO_PKEPLER))return TINYSECULAR_RANGE;
 n=sqrt(m->mu_m3_s2/a)/a;p=a*(1-ecc*ecc);c=cos(e->inclination_rad);
 factor=1.5*m->j2*n*pow(m->equatorial_radius_m/p,2);
 r.raan_rad_s=-factor*c;r.arg_periapsis_rad_s=.5*factor*(5*c*c-1);
 r.mean_anomaly_rad_s=n;
 if(profile==TINYSECULAR_AVERAGED_J2)r.mean_anomaly_rad_s+=.5*factor*sqrt(1-ecc*ecc)*(3*c*c-1);
 if(!isfinite(n)||n<=0||!isfinite(r.raan_rad_s)||!isfinite(r.arg_periapsis_rad_s)||!isfinite(r.mean_anomaly_rad_s))return TINYSECULAR_RANGE;
 *out=r;return TINYSECULAR_OK;
}
tinysecular_status tinysecular_propagate(const tinysecular_model *m,const tinysecular_elements *initial,tinysecular_profile profile,double elapsed,tinysecular_elements *out) {
 tinysecular_rates rates;tinysecular_elements r;double d0,d1,d2;tinysecular_status status;
 if(!out||!isfinite(elapsed))return TINYSECULAR_INVALID;
 status=tinysecular_j2_rates(m,initial,profile,&rates);if(status!=TINYSECULAR_OK)return status;
 d0=rates.raan_rad_s*elapsed;d1=rates.arg_periapsis_rad_s*elapsed;d2=rates.mean_anomaly_rad_s*elapsed;
 if(!isfinite(d0)||!isfinite(d1)||!isfinite(d2)||fabs(d0)>1e12||fabs(d1)>1e12||fabs(d2)>1e12)return TINYSECULAR_RANGE;
 r=*initial;r.raan_rad=wrap(wrap(r.raan_rad)+d0);r.arg_periapsis_rad=wrap(wrap(r.arg_periapsis_rad)+d1);r.mean_anomaly_rad=wrap(wrap(r.mean_anomaly_rad)+d2);
 *out=r;return TINYSECULAR_OK;
}
