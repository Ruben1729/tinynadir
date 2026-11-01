/* SPDX-License-Identifier: MIT */
#include "tinyomm/tinyomm.h"
#include "tinytime/tinytime.h"
#include <math.h>
#include <stdint.h>
#include <string.h>
static int space(char c){return c==' '||c=='\t'||c=='\r';}
static void trim(const char **s,size_t *n){while(*n&&space(**s)){++*s;--*n;}while(*n&&space((*s)[*n-1]))--*n;if(*n>=2&&((**s=='"'&&(*s)[*n-1]=='"')||(**s=='\''&&(*s)[*n-1]=='\''))){++*s;*n-=2;}}
static int equal(const char *s,size_t n,const char *text){return strlen(text)==n&&!memcmp(s,text,n);}
static int number(const char *s,size_t n,double *out){size_t i=0;uint64_t mantissa=0;unsigned kept=0;int sign=1,point=0,seen=0,fraction=0,dropped=0,exponent=0,esign=1;if(!n||n>128)return 0;if(s[i]=='+'||s[i]=='-'){if(s[i]=='-')sign=-1;++i;}for(;i<n;++i){char c=s[i];if(c=='.'&&!point){point=1;continue;}if(c=='e'||c=='E')break;if(c<'0'||c>'9')return 0;seen=1;if(point)++fraction;if(kept<19){mantissa=mantissa*10+(unsigned)(c-'0');++kept;}else ++dropped;}
    if(!seen)return 0;if(i<n){++i;if(i<n&&(s[i]=='+'||s[i]=='-')){if(s[i]=='-')esign=-1;++i;}if(i==n)return 0;for(;i<n;++i){if(s[i]<'0'||s[i]>'9'||exponent>1000)return 0;exponent=exponent*10+s[i]-'0';}}
    exponent=esign*exponent-fraction+dropped;if(exponent<-308||exponent>308)return 0;*out=(double)sign*(double)mantissa*pow(10,exponent);return isfinite(*out);
}
static int digits(const char *s,size_t n,unsigned *out){unsigned r=0;if(!n||n>9)return 0;for(size_t i=0;i<n;++i){if(s[i]<'0'||s[i]>'9')return 0;r=r*10+(unsigned)(s[i]-'0');}*out=r;return 1;}
static int epoch(const char *s,size_t n,double *first,double *fraction){unsigned year,month,day,hour,minute;double second;tinytime_jd jd;if(n&&s[n-1]=='Z')--n;if(n<19||n>48||s[4]!='-'||s[7]!='-'||s[10]!='T'||s[13]!=':'||s[16]!=':'||!digits(s,4,&year)||!digits(s+5,2,&month)||!digits(s+8,2,&day)||!digits(s+11,2,&hour)||!digits(s+14,2,&minute)||!number(s+17,n-17,&second))return 0;if(s[17]<'0'||s[17]>'9'||s[18]<'0'||s[18]>'9'||(n>19&&(s[19]!='.'||n==20)))return 0;for(size_t i=20;i<n;++i)if(s[i]<'0'||s[i]>'9')return 0;if(tinytime_calendar_to_jd((int)year,(int)month,(int)day,(int)hour,(int)minute,second,&jd)!=0)return 0;*first=jd.first;*fraction=jd.fraction;return 1;}
static int known(const char *s,size_t n){static const char *const keys[]={"EPOCH","MEAN_MOTION","ECCENTRICITY","INCLINATION","RA_OF_ASC_NODE","ARG_OF_PERICENTER","MEAN_ANOMALY","BSTAR","CENTER_NAME","REF_FRAME","TIME_SYSTEM","MEAN_ELEMENT_THEORY","EPHEMERIS_TYPE","NORAD_CAT_ID","CCSDS_OMM_VERS"};for(int i=0;i<15;++i)if(equal(s,n,keys[i]))return i;return -1;}
static int numeric_value(const char *s,size_t n,int key,double *out){size_t end=n;const char *unit=NULL;size_t unit_length=0;for(size_t i=0;i<n;++i)if(s[i]=='['){end=i;unit=s+i+1;if(n<i+2||s[n-1]!=']')return 0;unit_length=n-i-2;break;}while(end&&space(s[end-1]))--end;
    if(unit){const char *expected=key==1?"rev/day":((key>=3&&key<=6)?"deg":(key==7?"1/ER":"1"));if(!equal(unit,unit_length,expected))return 0;}
    return number(s,end,out);
}
tinyomm_status tinyomm_parse_fields(const tinyomm_field *fields,size_t count,int defaults,tinytle_elements *out){
    tinytle_elements g={0};unsigned mask=0;double *targets[8]={NULL,&g.mean_motion_rev_day,&g.eccentricity,&g.inclination_deg,&g.raan_deg,&g.arg_perigee_deg,&g.mean_anomaly_deg,&g.bstar};
    if(!out||!fields||!count||count>64||(defaults!=0&&defaults!=1))return TINYOMM_INVALID;
    for(size_t i=0;i<count;++i){const char *key=fields[i].name,*value=fields[i].value;size_t nk=fields[i].name_length,nv=fields[i].value_length;int k;unsigned integer;
        if(!key||!value||nk>64||nv>128)return TINYOMM_INVALID;trim(&key,&nk);trim(&value,&nv);k=known(key,nk);if(k<0)continue;if(mask&(1u<<k))return TINYOMM_INVALID;mask|=1u<<k;
        if(k==0){if(!epoch(value,nv,&g.epoch_jd1,&g.epoch_jd2))return TINYOMM_INVALID;}
        else if(k<8){if(!numeric_value(value,nv,k,targets[k]))return TINYOMM_INVALID;}
        else if(k<12){const char *expected=k==8?"EARTH":(k==9?"TEME":(k==10?"UTC":"SGP4"));if(!equal(value,nv,expected))return TINYOMM_UNSUPPORTED;}
        else if(k==12){if(!digits(value,nv,&integer))return TINYOMM_INVALID;if(integer!=0)return TINYOMM_UNSUPPORTED;g.ephemeris_type=0;}
        else if(k==13){if(!digits(value,nv,&integer))return TINYOMM_INVALID;g.catalog_id=integer;}
        else if(!equal(value,nv,"2.0")&&!equal(value,nv,"3.0"))return TINYOMM_UNSUPPORTED;
    }
    if((mask&255u)!=255u||(!defaults&&(mask&0xf00u)!=0xf00u))return TINYOMM_INVALID;
    if(g.inclination_deg<0||g.inclination_deg>180||g.raan_deg<0||g.raan_deg>=360||g.eccentricity<0||g.eccentricity>=1||g.arg_perigee_deg<0||g.arg_perigee_deg>=360||g.mean_anomaly_deg<0||g.mean_anomaly_deg>=360||g.mean_motion_rev_day<=0)return TINYOMM_INVALID;
    *out=g;return TINYOMM_OK;
}
tinyomm_status tinyomm_parse_kvn(const char *text,size_t length,int defaults,tinytle_elements *out){
    tinyomm_field fields[64];size_t count=0,pos=0;if(!text||!out||!length||length>8192)return TINYOMM_INVALID;
    while(pos<length){size_t start=pos,end,equal_pos;while(pos<length&&text[pos]!='\n')++pos;end=pos;if(pos<length)++pos;if(end-start>256)return TINYOMM_INVALID;equal_pos=start;while(equal_pos<end&&text[equal_pos]!='=')++equal_pos;
        if(equal_pos==end){const char *line=text+start;size_t n=end-start;trim(&line,&n);if(n&& !equal(line,n,"META_START")&&!equal(line,n,"META_STOP")&&!equal(line,n,"DATA_START")&&!equal(line,n,"DATA_STOP")&&!(n>=7&&!memcmp(line,"COMMENT",7)))return TINYOMM_INVALID;continue;}
        if(count==64)return TINYOMM_INVALID;fields[count].name=text+start;fields[count].name_length=equal_pos-start;fields[count].value=text+equal_pos+1;fields[count].value_length=end-equal_pos-1;++count;
    }
    return tinyomm_parse_fields(fields,count,defaults,out);
}
