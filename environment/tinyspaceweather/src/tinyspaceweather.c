/* SPDX-License-Identifier: MIT */
#include "tinyspaceweather/tinyspaceweather.h"
#include <math.h>
static const double ap_nodes[28]={0,2,3,4,5,6,7,9,12,15,18,22,27,32,39,48,56,67,80,94,111,132,154,179,207,236,300,400};
tinyspaceweather_status tinyspaceweather_kp_to_ap(double kp,double *ap){size_t i;double fraction,value;if(!ap||!isfinite(kp))return TINYSPACEWEATHER_INVALID;if(kp<0||kp>9)return TINYSPACEWEATHER_RANGE;i=(size_t)floor(kp*3);if(i==27){*ap=400;return TINYSPACEWEATHER_OK;}fraction=kp*3-(double)i;value=(1-fraction)*ap_nodes[i]+fraction*ap_nodes[i+1];*ap=value;return TINYSPACEWEATHER_OK;}
tinyspaceweather_status tinyspaceweather_ap_to_kp(double ap,double *kp){size_t i=0;if(!kp||!isfinite(ap))return TINYSPACEWEATHER_INVALID;if(ap<0||ap>400)return TINYSPACEWEATHER_RANGE;while(i<26&&ap>ap_nodes[i+1])++i;*kp=((double)i+(ap-ap_nodes[i])/(ap_nodes[i+1]-ap_nodes[i]))/3;return TINYSPACEWEATHER_OK;}
tinyspaceweather_status tinyspaceweather_interpolate(const tinyspaceweather_sample *table,size_t n,double t,tinyspaceweather_sample *out){
    size_t lo=0,hi;double fraction,span;tinyspaceweather_sample r;
    if(!table||!out||n<2||n>4096||!isfinite(t))return TINYSPACEWEATHER_INVALID;
    for(size_t i=0;i<n;++i){const tinyspaceweather_sample *s=&table[i];if(!isfinite(s->time_s)||!isfinite(s->f107_sfu)||s->f107_sfu<0||!isfinite(s->f107_average_sfu)||s->f107_average_sfu<0||!isfinite(s->ap)||s->ap<0||(i&&s->time_s<=table[i-1].time_s))return TINYSPACEWEATHER_INVALID;}
    if(t<table[0].time_s||t>table[n-1].time_s)return TINYSPACEWEATHER_RANGE;
    hi=n-1;while(hi-lo>1){size_t m=lo+(hi-lo)/2;if(t<table[m].time_s)hi=m;else lo=m;}
    span=table[hi].time_s-table[lo].time_s;fraction=(t-table[lo].time_s)/span;r.time_s=t;
    r.f107_sfu=(1-fraction)*table[lo].f107_sfu+fraction*table[hi].f107_sfu;
    r.f107_average_sfu=(1-fraction)*table[lo].f107_average_sfu+fraction*table[hi].f107_average_sfu;
    r.ap=(1-fraction)*table[lo].ap+fraction*table[hi].ap;
    if(!isfinite(span)||!isfinite(r.f107_sfu)||!isfinite(r.f107_average_sfu)||!isfinite(r.ap))return TINYSPACEWEATHER_ARITHMETIC;
    *out=r;return TINYSPACEWEATHER_OK;
}
