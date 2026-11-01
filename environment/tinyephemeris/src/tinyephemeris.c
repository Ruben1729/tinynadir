/* SPDX-License-Identifier: MIT */
#include "tinyephemeris/tinyephemeris.h"
#include "tinyinterpolation/tinyinterpolation.h"
#include <math.h>
#include <string.h>
tinyephemeris_status tinyephemeris_evaluate(const tinyephemeris_sample *samples,size_t n,double t,double position[3],double velocity[3]){
    size_t lo=0,hi;double p[3],v[3];
    if(!samples||!position||!velocity||position==velocity||n<2||n>4096||!isfinite(t))return TINYEPHEMERIS_INVALID;
    for(size_t i=0;i<n;++i){if(!isfinite(samples[i].time_s)||(i&&samples[i].time_s<=samples[i-1].time_s))return TINYEPHEMERIS_INVALID;for(size_t j=0;j<3;++j)if(!isfinite(samples[i].position_m[j])||!isfinite(samples[i].velocity_m_s[j]))return TINYEPHEMERIS_INVALID;}
    if(t<samples[0].time_s||t>samples[n-1].time_s)return TINYEPHEMERIS_RANGE;
    hi=n-1;while(hi-lo>1){size_t m=lo+(hi-lo)/2;if(t<samples[m].time_s)hi=m;else lo=m;}
    for(size_t j=0;j<3;++j)if(tinyinterpolation_hermite(samples[lo].time_s,samples[hi].time_s,samples[lo].position_m[j],samples[hi].position_m[j],samples[lo].velocity_m_s[j],samples[hi].velocity_m_s[j],t,&p[j],&v[j])!=TINYINTERPOLATION_OK)return TINYEPHEMERIS_ARITHMETIC;
    memcpy(position,p,sizeof p);memcpy(velocity,v,sizeof v);return TINYEPHEMERIS_OK;
}
