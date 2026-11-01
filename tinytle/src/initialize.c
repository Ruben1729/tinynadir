/* SPDX-License-Identifier: MIT */
#include "tinytle/tinytle.h"
#include <math.h>
static int valid(const tinytle_elements *g){
    const double values[]={g->epoch_jd1,g->epoch_jd2,g->inclination_deg,g->raan_deg,g->eccentricity,g->arg_perigee_deg,g->mean_anomaly_deg,g->mean_motion_rev_day,g->bstar};
    for(size_t i=0;i<9;++i)if(!isfinite(values[i]))return 0;
    return g->catalog_id<=999999999&&g->epoch_jd1>=2305447.5&&g->epoch_jd1<=5373483.5&&floor(g->epoch_jd1+.5)==g->epoch_jd1+.5&&g->epoch_jd2>=0&&g->epoch_jd2<1&&g->inclination_deg>=0&&g->inclination_deg<=180&&g->raan_deg>=0&&g->raan_deg<360&&g->eccentricity>=0&&g->eccentricity<1&&g->arg_perigee_deg>=0&&g->arg_perigee_deg<360&&g->mean_anomaly_deg>=0&&g->mean_anomaly_deg<360&&g->mean_motion_rev_day>0;
}
tinytle_status tinytle_initialize_elements(const tinytle_elements *g,tgp4_satellite *sat){
    tgp4_record r={0};const double deg2rad=TGP4_PI/180.0,rev_per_rad=1440.0/(2.0*TGP4_PI);
    if(!g||!sat||!valid(g))return TINYTLE_INVALID;if(g->ephemeris_type!=0)return TINYTLE_UNSUPPORTED;
    r.whichconst=TGP4_WGS72;r.jdsatepoch=g->epoch_jd1;r.jdsatepochF=g->epoch_jd2;r.bstar=g->bstar;r.ecco=g->eccentricity;r.argpo=g->arg_perigee_deg*deg2rad;r.inclo=g->inclination_deg*deg2rad;r.mo=g->mean_anomaly_deg*deg2rad;r.no_kozai=g->mean_motion_rev_day/rev_per_rad;r.nodeo=g->raan_deg*deg2rad;
    if(!tgp4_core_sgp4init('i',&r)||r.error||!isfinite(r.am)||!isfinite(r.nm))return TINYTLE_PROPAGATION_ERROR;
    sat->record=r;sat->initialized=1;return TINYTLE_OK;
}
tinytle_status tinytle_initialize(const char *a,size_t na,const char *b,size_t nb,tgp4_satellite *sat){tinytle_elements g;tinytle_status status=tinytle_parse(a,na,b,nb,&g);if(status!=TINYTLE_OK)return status;return tinytle_initialize_elements(&g,sat);}
