/* SPDX-License-Identifier: MIT */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL %s:%d %s\n",__FILE__,__LINE__,#c); return 1; } } while (0)

#include "tinyomm/tinyomm.h"
#include "data/reference.h"
#include "data/decimal.h"
#include "data/exponents.h"
static tinyomm_status parse_bstar(const char *text,tinytle_elements *g){
 const char *names[8]={"EPOCH","MEAN_MOTION","ECCENTRICITY","INCLINATION","RA_OF_ASC_NODE","ARG_OF_PERICENTER","MEAN_ANOMALY","BSTAR"};
 const char *values[8]={"2026-11-01T12:00:00Z","15",".001","51","10","20","30",text};tinyomm_field fields[8];
 for(size_t i=0;i<8;++i){fields[i].name=names[i];fields[i].name_length=strlen(names[i]);fields[i].value=values[i];fields[i].value_length=strlen(values[i]);}
 return tinyomm_parse_fields(fields,8,1,g);
}
int main(void){double maxp=0,maxv=0;
 for(size_t i=0;i<sizeof cases/sizeof cases[0];++i){const reference_case *c=&cases[i];tinytle_elements g;tgp4_satellite sat={0};tgp4_state state;CHECK(tinyomm_parse_kvn(c->kvn,strlen(c->kvn),0,&g)==0);CHECK(g.catalog_id==c->catalog_id);CHECK(tinyomm_initialize_kvn(c->kvn,strlen(c->kvn),0,&sat)==0);for(size_t j=0;j<3;++j){CHECK(tgp4_propagate(&sat,c->time[j],&state)==0);for(size_t k=0;k<3;++k){double ep=fabs(state.position_km[k]-c->states[j][k])*1000,ev=fabs(state.velocity_km_s[k]-c->states[j][3+k])*1000;maxp=fmax(maxp,ep);maxv=fmax(maxv,ev);CHECK(ep<.001);CHECK(ev<1e-5);}}}
 {char buffer[2048];tinytle_elements g={0},saved=g;size_t n=strlen(cases[0].kvn);CHECK(n+80<sizeof buffer);memcpy(buffer,cases[0].kvn,n);memcpy(buffer+n,"EPHEMERIS_TYPE = 4\n",20);CHECK(tinyomm_parse_kvn(buffer,n+19,0,&g)!=0);CHECK(!memcmp(&g,&saved,sizeof g));
 const char *minimal="EPOCH = 2026-11-01T12:00:00Z\nMEAN_MOTION = 15\nECCENTRICITY = .001\nINCLINATION = 51\nRA_OF_ASC_NODE = 10\nARG_OF_PERICENTER = 20\nMEAN_ANOMALY = 30\nBSTAR = 1e-5\n";CHECK(tinyomm_parse_kvn(minimal,strlen(minimal),0,&g)==TINYOMM_INVALID);CHECK(tinyomm_parse_kvn(minimal,strlen(minimal),1,&g)==0);CHECK(tinyomm_parse_kvn(NULL,10,1,&g)==TINYOMM_INVALID);
 tinyomm_field fields[]={{"EPOCH",5,"2000-01-01T12:00:00",19},{"MEAN_MOTION",11,"15",2}};CHECK(tinyomm_parse_fields(fields,2,1,&g)==TINYOMM_INVALID);}
 for(size_t i=0;i<sizeof decimals/sizeof decimals[0];++i){
  const char *names[8]={"EPOCH","MEAN_MOTION","ECCENTRICITY","INCLINATION","RA_OF_ASC_NODE","ARG_OF_PERICENTER","MEAN_ANOMALY","BSTAR"};
  const char *values[8]={"2026-11-01T12:00:00Z","15",".001","51","10","20","30",decimals[i].text};
  tinyomm_field fields[8];tinytle_elements g;
  for(size_t j=0;j<8;++j){fields[j].name=names[j];fields[j].name_length=strlen(names[j]);fields[j].value=values[j];fields[j].value_length=strlen(values[j]);}
  CHECK(tinyomm_parse_fields(fields,8,1,&g)==0);CHECK(fabs(g.bstar-decimals[i].expected)<=4e-15*fmax(1e-300,fabs(decimals[i].expected)));
 }
 for(size_t i=0;i<sizeof exponent_cases/sizeof exponent_cases[0];++i){
  tinytle_elements g;double expected=exponent_cases[i].expected;
  double ulp=expected==0?nextafter(0,1):fabs(expected-nextafter(expected,0));
  CHECK(parse_bstar(exponent_cases[i].text,&g)==TINYOMM_OK);
  CHECK(isfinite(g.bstar)&&fabs(g.bstar-expected)<=3*ulp);
  if(expected!=0)CHECK(g.bstar!=0);
  if(expected==0){CHECK(g.bstar==0);CHECK(!!signbit(g.bstar)==!!signbit(expected));}
 }
 {static const char *const invalid[]={"1e309","-1e309","1.7976931348623159e308","1.7976931348623158080e308","1.797693134862315807938e308","1e","1e+","nan","inf","1e--5","1.2.3","e5"};tinytle_elements g,saved;
  CHECK(parse_bstar("1e-308",&g)==TINYOMM_OK);double equivalent=g.bstar;
  CHECK(parse_bstar("1.0e-308",&g)==TINYOMM_OK&&g.bstar==equivalent);
  CHECK(parse_bstar("1000000000000000000e-326",&g)==TINYOMM_OK&&g.bstar==equivalent);
  memcpy(&saved,&g,sizeof g);
  for(size_t i=0;i<sizeof invalid/sizeof invalid[0];++i){CHECK(parse_bstar(invalid[i],&g)==TINYOMM_INVALID);CHECK(!memcmp(&saved,&g,sizeof g));}
 }
 {const char *extreme="EPOCH=2026-11-01T12:00:00Z\nMEAN_MOTION=1e300\nECCENTRICITY=.001\nINCLINATION=51\nRA_OF_ASC_NODE=10\nARG_OF_PERICENTER=20\nMEAN_ANOMALY=30\nBSTAR=0\n";
  tgp4_satellite sat={0},saved;CHECK(tinyomm_initialize_kvn(cases[0].kvn,strlen(cases[0].kvn),0,&sat)==TINYOMM_OK);memcpy(&saved,&sat,sizeof sat);
  CHECK(tinyomm_initialize_kvn(extreme,strlen(extreme),1,&sat)==TINYOMM_PROPAGATION_ERROR);CHECK(!memcmp(&sat,&saved,sizeof sat));
 }
 printf("tinyomm: independent OMM-to-SGP4 max %.12g m %.12g m/s\n",maxp,maxv);return 0;}
