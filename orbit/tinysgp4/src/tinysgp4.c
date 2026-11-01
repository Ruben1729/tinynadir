#include "tinysgp4/tinysgp4.h"
#include <stdbool.h>
#include <stdint.h>
#include <float.h>
#include "tinysgp4/math_compat.h"

/* C99 compile-time assertion: an unsupported double format gives a negative size. */
typedef char tgp4_requires_binary64_double[
    (sizeof(double) == 8 && DBL_MANT_DIG == 53 && DBL_MAX_EXP == 1024) ? 1 : -1];

static bool line_valid(const char* p,size_t n,char label){
  if(!p || n<69 || p[0]!=label || p[1]!=' ')return false;
  for(size_t i=69;i<n;++i)if(p[i]!='\r'&&p[i]!='\n')return false;
  int sum=0;
  for(int i=0;i<68;++i){if(p[i]>='0'&&p[i]<='9')sum+=p[i]-'0';else if(p[i]=='-')++sum;}
  return p[68]>='0' && p[68]<='9' && sum%10==p[68]-'0';
}
static bool decimal(const char* p,size_t n,double* out,bool implied){
  int sign=1;size_t i=0;while(i<n&&p[i]==' ')++i;
  if(i<n&&(p[i]=='+'||p[i]=='-')){if(p[i]=='-')sign=-1;++i;}
  int64_t value=0,den=1;bool point=implied,seen=false;
  for(;i<n;++i){char c=p[i];
    if(c==' '){while(i<n&&p[i]==' ')++i;if(i!=n)return false;break;}
    if(c=='.'&&!point){point=true;continue;}
    if(c<'0'||c>'9')return false;
    if(value>INT64_MAX/10 || den>INT64_MAX/10)return false;
    value=value*10+c-'0';if(point)den*=10;seen=true;
  }
  if(!seen)return false;
  *out=(double)sign*(double)value/(double)den;return TGP4_ISFINITE(*out);
}
static bool exponential(const char* p,double* out){
  bool blank=true;for(int i=0;i<8;++i)if(p[i]!=' ')blank=false;
  if(blank){*out=0;return true;}
  if((p[0]!=' '&&p[0]!='+'&&p[0]!='-')||(p[6]!='+'&&p[6]!='-')||p[7]<'0'||p[7]>'9')return false;
  int64_t v=0;for(int i=1;i<=5;++i){if(p[i]<'0'||p[i]>'9')return false;v=v*10+p[i]-'0';}
  int e=(p[7]-'0')*(p[6]=='-'?-1:1);
  double mantissa=(double)v/100000.0;
  if(p[0]=='-')mantissa=-mantissa;
  *out=mantissa*pow(10.0,e);return TGP4_ISFINITE(*out);
}
static int days_before_year(int y){--y;return 365*y+y/4-y/100+y/400;}
tgp4_status tgp4_initialize(tgp4_satellite* sat,const char* a,size_t na,const char* b,size_t nb){
 if(!sat || !line_valid(a,na,'1')||!line_valid(b,nb,'2'))return TGP4_INVALID_TLE;
 for(int k=2;k<7;++k)if(a[k]!=b[k] || a[k]<'0'||a[k]>'9')return TGP4_INVALID_TLE;
 if(a[18]<'0' || a[18]>'9' || a[19]<'0' || a[19]>'9')return TGP4_INVALID_TLE;
 if(a[62]!='0')return TGP4_UNSUPPORTED_EPHEMERIS;
 double yy,day,ecc,inc,node,arg,mean,motion,bstar;
 if(!decimal(a+18,2,&yy,false)||!decimal(a+20,12,&day,false)||!exponential(a+53,&bstar)||
 !decimal(b+8,8,&inc,false)||!decimal(b+17,8,&node,false)||!decimal(b+26,7,&ecc,true)||
 !decimal(b+34,8,&arg,false)||!decimal(b+43,8,&mean,false)||!decimal(b+52,11,&motion,false))return TGP4_INVALID_TLE;
 int year=(int)yy;year+=year<57?2000:1900;
 int ndays=days_before_year(year+1)-days_before_year(year);
 if(day<1 || day>=ndays+1 || inc<0 || inc>180 || node<0 || node>=360 ||
 arg<0 || arg>=360 || mean<0 || mean>=360 || ecc<0 || ecc>=1 || motion<=0)return TGP4_INVALID_TLE;
 tgp4_record r={0};r.whichconst=TGP4_WGS72;
 double whole=floor(day);
 r.jdsatepoch=2433281.5+(days_before_year(year)-days_before_year(1950))+whole;
 r.jdsatepochF=day-whole;r.bstar=bstar;r.ecco=ecc;
 /* Preserve Vallado twoline2rv's operation order. Algebraically equivalent
    conversions can diverge in ill-conditioned/high-drag verification cases. */
 const double deg2rad=TGP4_PI/180.0;
 const double rev_per_rad=1440.0/(2.0*TGP4_PI);
 r.argpo=arg*deg2rad;r.inclo=inc*deg2rad;r.mo=mean*deg2rad;
 r.no_kozai=motion/rev_per_rad;r.nodeo=node*deg2rad;
 bool ok=tgp4_core_sgp4init('i',&r);
 if(!ok || r.error)return TGP4_PROPAGATION_ERROR;
 if(!TGP4_ISFINITE(r.am) || !TGP4_ISFINITE(r.nm))return TGP4_ARITHMETIC_ERROR;
 sat->record=r;sat->initialized=1;return TGP4_OK;
}
tgp4_status tgp4_propagate(tgp4_satellite* sat,double minutes,tgp4_state* output){
 if(!sat || !output || !sat->initialized)return TGP4_UNINITIALIZED;
 if(!TGP4_ISFINITE(minutes) || minutes < -10080 || minutes>10080)return TGP4_TIME_OUT_OF_RANGE;
 tgp4_record next=sat->record;tgp4_state result;
 bool ok=tgp4_core_sgp4(&next,minutes,result.position_km,result.velocity_km_s);
 if(!ok || next.error)return TGP4_PROPAGATION_ERROR;
 for(int i=0;i<3;++i)
  if(!TGP4_ISFINITE(result.position_km[i]) || !TGP4_ISFINITE(result.velocity_km_s[i]))return TGP4_ARITHMETIC_ERROR;
 sat->record=next;*output=result;return TGP4_OK;
}
