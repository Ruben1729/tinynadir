/* SPDX-License-Identifier: MIT */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL %s:%d %s\n",__FILE__,__LINE__,#c); return 1; } } while (0)

#include "tinytle/tinytle.h"
#include "data/reference.h"
static void checksum(char *s){int sum=0;for(size_t i=0;i<68;++i){if(s[i]>='0'&&s[i]<='9')sum+=s[i]-'0';else if(s[i]=='-')++sum;}s[68]=(char)('0'+sum%10);}
int main(void){double maxp=0,maxv=0;
 for(size_t i=0;i<sizeof cases/sizeof cases[0];++i){const reference_case *c=&cases[i];tinytle_elements g;tgp4_satellite sat={0};tgp4_state state;CHECK(tinytle_parse(c->a,69,c->b,69,&g)==0);CHECK(g.catalog_id==c->id);CHECK(g.epoch_jd1==c->epoch[0]);CHECK(fabs(g.epoch_jd2-c->epoch[1])<5e-14);CHECK(tinytle_initialize_elements(&g,&sat)==0);for(size_t j=0;j<3;++j){CHECK(tgp4_propagate(&sat,c->time[j],&state)==0);for(size_t k=0;k<3;++k){double ep=fabs(state.position_km[k]-c->states[j][k])*1000,ev=fabs(state.velocity_km_s[k]-c->states[j][3+k])*1000;maxp=fmax(maxp,ep);maxv=fmax(maxv,ev);CHECK(ep<.001);CHECK(ev<1e-5);}}}
 {char a[70],b[70];tinytle_elements g={0},saved=g;tgp4_satellite sat={0};memcpy(a,cases[0].a,70);memcpy(b,cases[0].b,70);for(size_t n=0;n<69;++n)CHECK(tinytle_parse(a,n,b,69,&g)==TINYTLE_INVALID);a[68]=a[68]=='9'?'0':a[68]+1;CHECK(tinytle_parse(a,69,b,69,&g)==TINYTLE_INVALID);CHECK(!memcmp(&g,&saved,sizeof g));memcpy(a,cases[0].a,70);a[62]='4';checksum(a);CHECK(tinytle_parse(a,69,b,69,&g)==TINYTLE_UNSUPPORTED);memcpy(a,cases[0].a,70);memcpy(a+2,"I0000",5);memcpy(b+2,"I0000",5);checksum(a);checksum(b);CHECK(tinytle_parse(a,69,b,69,&g)==TINYTLE_INVALID);g.epoch_jd1=NAN;CHECK(tinytle_initialize_elements(&g,&sat)==TINYTLE_INVALID);}
 printf("tinytle: independent extended-TLE SGP4 max %.12g m %.12g m/s\n",maxp,maxv);return 0;}
