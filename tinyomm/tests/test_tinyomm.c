/* SPDX-License-Identifier: MIT */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL %s:%d %s\n",__FILE__,__LINE__,#c); return 1; } } while (0)

#include "tinyomm/tinyomm.h"
#include "data/reference.h"
int main(void){double maxp=0,maxv=0;
 for(size_t i=0;i<sizeof cases/sizeof cases[0];++i){const reference_case *c=&cases[i];tinytle_elements g;tgp4_satellite sat={0};tgp4_state state;CHECK(tinyomm_parse_kvn(c->kvn,strlen(c->kvn),0,&g)==0);CHECK(g.catalog_id==c->catalog_id);CHECK(tinyomm_initialize_kvn(c->kvn,strlen(c->kvn),0,&sat)==0);for(size_t j=0;j<3;++j){CHECK(tgp4_propagate(&sat,c->time[j],&state)==0);for(size_t k=0;k<3;++k){double ep=fabs(state.position_km[k]-c->states[j][k])*1000,ev=fabs(state.velocity_km_s[k]-c->states[j][3+k])*1000;maxp=fmax(maxp,ep);maxv=fmax(maxv,ev);CHECK(ep<.001);CHECK(ev<1e-5);}}}
 {char buffer[2048];tinytle_elements g={0},saved=g;size_t n=strlen(cases[0].kvn);CHECK(n+80<sizeof buffer);memcpy(buffer,cases[0].kvn,n);memcpy(buffer+n,"EPHEMERIS_TYPE = 4\n",20);CHECK(tinyomm_parse_kvn(buffer,n+19,0,&g)!=0);CHECK(!memcmp(&g,&saved,sizeof g));
 const char *minimal="EPOCH = 2026-11-01T12:00:00Z\nMEAN_MOTION = 15\nECCENTRICITY = .001\nINCLINATION = 51\nRA_OF_ASC_NODE = 10\nARG_OF_PERICENTER = 20\nMEAN_ANOMALY = 30\nBSTAR = 1e-5\n";CHECK(tinyomm_parse_kvn(minimal,strlen(minimal),0,&g)==TINYOMM_INVALID);CHECK(tinyomm_parse_kvn(minimal,strlen(minimal),1,&g)==0);CHECK(tinyomm_parse_kvn(NULL,10,1,&g)==TINYOMM_INVALID);
 tinyomm_field fields[]={{"EPOCH",5,"2000-01-01T12:00:00",19},{"MEAN_MOTION",11,"15",2}};CHECK(tinyomm_parse_fields(fields,2,1,&g)==TINYOMM_INVALID);}
 printf("tinyomm: independent OMM-to-SGP4 max %.12g m %.12g m/s\n",maxp,maxv);return 0;}
