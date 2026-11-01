/* SPDX-License-Identifier: MIT */
#include "tinytle/tinytle.h"
#include "tinytime/tinytime.h"
#include <math.h>
#include <stdint.h>
static int line_valid(const char *s,size_t n,char tag){int sum=0;if(!s||n<69||n>71||s[0]!=tag||s[1]!=' ')return 0;for(size_t i=69;i<n;++i)if(s[i]!='\r'&&s[i]!='\n')return 0;for(size_t i=0;i<68;++i){if(s[i]>='0'&&s[i]<='9')sum+=s[i]-'0';else if(s[i]=='-')++sum;}return s[68]>='0'&&s[68]<='9'&&sum%10==s[68]-'0';}
static int decimal(const char *s,size_t n,int implied,double *out){size_t i=0;uint64_t value=0,den=1;int sign=1,point=implied,seen=0;while(i<n&&s[i]==' ')++i;if(i<n&&(s[i]=='+'||s[i]=='-')){if(s[i]=='-')sign=-1;++i;}for(;i<n;++i){unsigned digit;if(s[i]==' '){while(i<n&&s[i]==' ')++i;if(i!=n)return 0;break;}if(s[i]=='.'&&!point){point=1;continue;}if(s[i]<'0'||s[i]>'9')return 0;digit=(unsigned)(s[i]-'0');if(value>(UINT64_MAX-digit)/10||(point&&den>UINT64_MAX/10))return 0;value=value*10+digit;if(point)den*=10;seen=1;}if(!seen)return 0;*out=(double)sign*(double)value/(double)den;return isfinite(*out);}
static int exponential(const char *s,double *out){double mantissa;int exponent;size_t blanks=0;while(blanks<8&&s[blanks]==' ')++blanks;if(blanks==8){*out=0;return 1;}if((s[0]!=' '&&s[0]!='+'&&s[0]!='-')||(s[6]!='+'&&s[6]!='-')||s[7]<'0'||s[7]>'9')return 0;for(size_t i=1;i<=5;++i)if(s[i]<'0'||s[i]>'9')return 0;if(!decimal(s,6,1,&mantissa))return 0;exponent=(s[7]-'0')*(s[6]=='-'?-1:1);*out=mantissa*pow(10,exponent);return isfinite(*out);}
static int identifier(const char *s,uint32_t *id){
    static const char alphabet[]="ABCDEFGHJKLMNPQRSTUVWXYZ";
    uint32_t value=0;size_t i=0;
    if(s[0]==' '||(s[0]>='0'&&s[0]<='9')){
        while(i<5&&s[i]==' ')++i;
        if(i==5)return 0;
        for(;i<5;++i){if(s[i]<'0'||s[i]>'9')return 0;value=value*10+(unsigned)(s[i]-'0');}
    }else{
        int found=0;uint32_t suffix=0;
        for(unsigned j=0;j<24;++j)if(s[0]==alphabet[j]){value=(j+10)*10000;found=1;break;}
        if(!found)return 0;
        for(i=1;i<5;++i){if(s[i]<'0'||s[i]>'9')return 0;suffix=suffix*10+(unsigned)(s[i]-'0');}
        value+=suffix;
    }
    *id=value;return 1;
}
tinytle_status tinytle_parse(const char *a,size_t na,const char *b,size_t nb,tinytle_elements *out){
    tinytle_elements g;uint32_t id2;double yy,day;int year;double maximum;tinytime_jd jan,next,epoch;
    if(!out||!line_valid(a,na,'1')||!line_valid(b,nb,'2')||!identifier(a+2,&g.catalog_id)||!identifier(b+2,&id2)||id2!=g.catalog_id)return TINYTLE_INVALID;
    if(a[62]!='0')return TINYTLE_UNSUPPORTED;
    g.ephemeris_type=0;
    if(a[18]<'0'||a[18]>'9'||a[19]<'0'||a[19]>'9'||!decimal(a+18,2,0,&yy)||!decimal(a+20,12,0,&day)||!exponential(a+53,&g.bstar)||!decimal(b+8,8,0,&g.inclination_deg)||!decimal(b+17,8,0,&g.raan_deg)||!decimal(b+26,7,1,&g.eccentricity)||!decimal(b+34,8,0,&g.arg_perigee_deg)||!decimal(b+43,8,0,&g.mean_anomaly_deg)||!decimal(b+52,11,0,&g.mean_motion_rev_day))return TINYTLE_INVALID;
    year=(int)yy;year+=year<57?2000:1900;
    if(tinytime_calendar_to_jd(year,1,1,0,0,0,&jan)!=0||tinytime_calendar_to_jd(year+1,1,1,0,0,0,&next)!=0)return TINYTLE_INVALID;
    maximum=next.first-jan.first;
    if(day<0||day>=maximum+1||g.inclination_deg<0||g.inclination_deg>180||g.raan_deg<0||g.raan_deg>=360||g.eccentricity<0||g.eccentricity>=1||g.arg_perigee_deg<0||g.arg_perigee_deg>=360||g.mean_anomaly_deg<0||g.mean_anomaly_deg>=360||g.mean_motion_rev_day<=0)return TINYTLE_INVALID;
    if(tinytime_add_seconds(jan,(floor(day)-1)*86400,&epoch)!=0)return TINYTLE_INVALID;
    g.epoch_jd1=epoch.first;g.epoch_jd2=day-floor(day);*out=g;return TINYTLE_OK;
}
