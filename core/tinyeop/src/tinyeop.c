/* SPDX-License-Identifier: MIT */
#include "tinyeop/tinyeop.h"
#include <math.h>

typedef double (*getter)(const void*,size_t,size_t);
static double eop_value(const void *data,size_t i,size_t column) {
    const tinyeop_record *r=data;return column?r[i].value[column-1]:r[i].mjd_utc;
}
static double xys_value(const void *data,size_t i,size_t column) {
    const tinyeop_xys_record *r=data;
    if(column==0)return r[i].mjd_tt;
    if(column==1)return r[i].x_rad;
    if(column==2)return r[i].y_rad;
    return r[i].s_rad;
}
static tinyeop_status window(const void *base,getter get,size_t n,double t,tinyeop_method method,size_t *start,size_t *width) {
    size_t i,low=0,high,k;
    if(!base||!isfinite(t))return TINYEOP_INVALID;
    if(n<2||n>TINYEOP_MAX_SAMPLES||(method!=TINYEOP_LINEAR&&method!=TINYEOP_CUBIC)||(method==TINYEOP_CUBIC&&n<4))return TINYEOP_RANGE;
    for(i=0;i<n;i++)if(!isfinite(get(base,i,0))||(i&&get(base,i,0)<=get(base,i-1,0)))return TINYEOP_TABLE;
    if(t<get(base,0,0)||t>get(base,n-1,0))return TINYEOP_RANGE;
    high=n-1;while(high-low>1){k=low+(high-low)/2;if(t<get(base,k,0))high=k;else low=k;}
    if(method==TINYEOP_LINEAR){*start=low;*width=2;}
    else {*start=low?low-1:0;if(*start>n-4)*start=n-4;*width=4;}
    return TINYEOP_OK;
}
static double interpolate(const void *base,getter get,size_t first,size_t width,size_t column,double t) {
    size_t i,j;double result=0;
    for(i=0;i<width;i++){
        double term=get(base,first+i,column);
        for(j=0;j<width;j++)if(j!=i)term*=(t-get(base,first+j,0))/(get(base,first+i,0)-get(base,first+j,0));
        result+=term;
    }
    return result;
}
tinyeop_status tinyeop_sample(const tinyeop_record *records,size_t n,double t,double dat,tinyeop_method method,tinyeop_result *out) {
    size_t first,width,i,j;double v[6];tinyeop_result r;tinyeop_status s;
    if(!records||!out||!isfinite(dat))return TINYEOP_INVALID;
    if(fabs(dat)>1000)return TINYEOP_RANGE;
    s=window(records,eop_value,n,t,method,&first,&width);if(s!=TINYEOP_OK)return s;
    for(i=0;i<n;i++)for(j=0;j<6;j++)if(!isfinite(records[i].value[j]))return TINYEOP_TABLE;
    for(j=0;j<6;j++){v[j]=interpolate(records,eop_value,first,width,j+1,t);if(!isfinite(v[j]))return TINYEOP_TABLE;}
    r.xp_rad=v[0];r.yp_rad=v[1];r.dut1_s=v[2]+dat;r.lod_s=v[3];r.dx_rad=v[4];r.dy_rad=v[5];
    *out=r;return TINYEOP_OK;
}
tinyeop_status tinyeop_sample_xys(const tinyeop_xys_record *records,size_t n,double t,tinyeop_method method,double out[3]) {
    size_t first,width,i,j;double v[3];tinyeop_status s;
    if(!records||!out)return TINYEOP_INVALID;
    s=window(records,xys_value,n,t,method,&first,&width);if(s!=TINYEOP_OK)return s;
    for(i=0;i<n;i++)if(!isfinite(records[i].x_rad)||!isfinite(records[i].y_rad)||!isfinite(records[i].s_rad))return TINYEOP_TABLE;
    for(j=0;j<3;j++){v[j]=interpolate(records,xys_value,first,width,j+1,t);if(!isfinite(v[j]))return TINYEOP_TABLE;}
    out[0]=v[0];out[1]=v[1];out[2]=v[2];return TINYEOP_OK;
}
