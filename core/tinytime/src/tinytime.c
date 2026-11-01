/* SPDX-License-Identifier: MIT */
#include "tinytime/tinytime.h"
#include <math.h>
#include <stddef.h>

static int leap(int y) { return y%4==0 && (y%100!=0 || y%400==0); }
static int month_days(int y,int m) {
    static const int d[]={31,28,31,30,31,30,31,31,30,31,30,31};
    return d[m-1]+(m==2 && leap(y));
}
tinytime_status tinytime_normalize(double a,double b,tinytime_jd *out) {
    tinytime_jd r; double whole;
    if(!out || !isfinite(a) || !isfinite(b))return TINYTIME_INVALID;
    if(fabs(a)>1e8 || fabs(b)>1e8)return TINYTIME_RANGE;
    r.first=floor(a-0.5)+0.5; r.fraction=(a-r.first)+b;
    whole=floor(r.fraction);r.first+=whole;r.fraction-=whole;
    if(fabs(r.first)>1e8)return TINYTIME_RANGE;
    *out=r;return TINYTIME_OK;
}
tinytime_status tinytime_calendar_to_jd(int y,int m,int d,int h,int min,double s,tinytime_jd *out) {
    long a,yy,mm,jdn; tinytime_jd r;
    if(!out || !isfinite(s))return TINYTIME_INVALID;
    if(y<1600 || y>9999 || m<1 || m>12 || d<1 || d>month_days(y,m) || h<0 || h>23 || min<0 || min>59 || s<0 || s>=60)return TINYTIME_RANGE;
    a=(14-m)/12; yy=y+4800-a;mm=m+12*a-3;
    jdn=d+(153*mm+2)/5+365*yy+yy/4-yy/100+yy/400-32045;
    r.first=(double)jdn-0.5;r.fraction=(h*3600.0+min*60.0+s)/86400.0;
    *out=r;return TINYTIME_OK;
}
tinytime_status tinytime_jd_to_calendar(tinytime_jd jd,int *y,int *m,int *d,int *h,int *min,double *s) {
    tinytime_jd r;long a,b,c,dd,e,mm;int yr,mo,day,hr,mn;double seconds;
    if(!y||!m||!d||!h||!min||!s)return TINYTIME_INVALID;
    if(tinytime_normalize(jd.first,jd.fraction,&r)!=TINYTIME_OK)return TINYTIME_INVALID;
    if(r.first<2305447.5 || r.first>=5373484.5)return TINYTIME_RANGE;
    a=(long)floor(r.first+0.5)+32044;b=(4*a+3)/146097;c=a-(146097*b)/4;
    dd=(4*c+3)/1461;e=c-1461*dd/4;mm=(5*e+2)/153;
    day=(int)(e-(153*mm+2)/5+1);mo=(int)(mm+3-12*(mm/10));yr=(int)(100*b+dd-4800+mm/10);
    seconds=r.fraction*86400.0;hr=(int)(seconds/3600.0);seconds-=hr*3600.0;mn=(int)(seconds/60.0);seconds-=mn*60.0;
    *y=yr;*m=mo;*d=day;*h=hr;*min=mn;*s=seconds;return TINYTIME_OK;
}
tinytime_status tinytime_add_seconds(tinytime_jd jd,double seconds,tinytime_jd *out) {
    if(!isfinite(seconds))return TINYTIME_INVALID;
    return tinytime_normalize(jd.first,jd.fraction+seconds/86400.0,out);
}
tinytime_status tinytime_difference_seconds(tinytime_jd a,tinytime_jd b,double *seconds) {
    double v;
    if(!seconds || !isfinite(a.first) || !isfinite(a.fraction) || !isfinite(b.first) || !isfinite(b.fraction))return TINYTIME_INVALID;
    v=((a.first-b.first)+(a.fraction-b.fraction))*86400.0;
    if(!isfinite(v))return TINYTIME_RANGE;
    *seconds=v;return TINYTIME_OK;
}
tinytime_status tinytime_decimal_year(tinytime_jd jd,double *year) {
    int y,m,d,h,min;double s,days;tinytime_jd start;
    if(!year)return TINYTIME_INVALID;
    if(tinytime_jd_to_calendar(jd,&y,&m,&d,&h,&min,&s)!=TINYTIME_OK)return TINYTIME_RANGE;
    if(tinytime_calendar_to_jd(y,1,1,0,0,0,&start)!=TINYTIME_OK)return TINYTIME_RANGE;
    if(tinytime_difference_seconds(jd,start,&days)!=TINYTIME_OK)return TINYTIME_RANGE;
    *year=y+days/(86400.0*(365+leap(y)));return TINYTIME_OK;
}
tinytime_status tinytime_utc_to_tai(tinytime_jd jd,double offset,tinytime_jd *out) {
    if(!isfinite(offset) || fabs(offset)>1000)return TINYTIME_RANGE;
    return tinytime_add_seconds(jd,offset,out);
}
tinytime_status tinytime_utc_to_ut1(tinytime_jd jd,double offset,tinytime_jd *out) {
    if(!isfinite(offset) || fabs(offset)>100)return TINYTIME_RANGE;
    return tinytime_add_seconds(jd,offset,out);
}
tinytime_status tinytime_tai_to_tt(tinytime_jd jd,tinytime_jd *out) {return tinytime_add_seconds(jd,32.184,out);}
tinytime_status tinytime_tt_to_tai(tinytime_jd jd,tinytime_jd *out) {return tinytime_add_seconds(jd,-32.184,out);}
tinytime_status tinytime_tt_to_tdb(tinytime_jd jd,tinytime_jd *out) {
    double t,g,offset;
    if(!isfinite(jd.first)||!isfinite(jd.fraction))return TINYTIME_INVALID;
    t=((jd.first-2451545.0)+jd.fraction)/36525.0;
    if(t < -1.0 || t > 1.0)return TINYTIME_RANGE;
    g=(357.5277233+35999.05034*t)*(3.14159265358979323846/180.0);
    offset=0.001657*sin(g)+0.000022*sin(2*g);
    return tinytime_add_seconds(jd,offset,out);
}
