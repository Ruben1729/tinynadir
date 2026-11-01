/* SPDX-License-Identifier: MIT */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL %s:%d %s\n",__FILE__,__LINE__,#c); return 1; } } while (0)

#include "tinycowell/tinycowell.h"
#include "data/reference.h"
static int drag(double t,const tinycowell_state *s,void *ctx,double a[3]){double k=*(double*)ctx;(void)t;for(size_t i=0;i<3;++i)a[i]=-k*s->velocity_m_s[i];return 0;}
static int fail(double t,const tinycowell_state *s,void *ctx,double a[3]){(void)t;(void)s;(void)ctx;(void)a;return 1;}
int main(void){double w[54],maxp=0,maxv=0;tinyrkf45_options o={1e-9,1e-12,1,1e-9,30,100000};tinyrkf45_stats stats;
    for(size_t i=0;i<sizeof cases/sizeof cases[0];++i){const reference_case *c=&cases[i];tinycowell_state initial,out;double coefficient=c->drag;const tinycowell_model model={c->mu,c->drag?drag:NULL,&coefficient};for(size_t j=0;j<3;++j){initial.position_m[j]=c->initial[j];initial.velocity_m_s[j]=c->initial[3+j];}
        CHECK(tinycowell_adaptive(&model,0,c->time,&initial,&o,&out,&stats,w,54)==0);
        for(size_t j=0;j<3;++j){maxp=fmax(maxp,fabs(out.position_m[j]-c->expected[j]));maxv=fmax(maxv,fabs(out.velocity_m_s[j]-c->expected[3+j]));CHECK(fabs(out.position_m[j]-c->expected[j])<.002);CHECK(fabs(out.velocity_m_s[j]-c->expected[3+j])<3e-6);}
        CHECK(tinycowell_fixed(&model,0,c->time,&initial,1,10000,&out,w,54)==0);for(size_t j=0;j<3;++j){CHECK(fabs(out.position_m[j]-c->expected[j])<.003);CHECK(fabs(out.velocity_m_s[j]-c->expected[3+j])<5e-6);}}
    {tinycowell_state s={{7000000,0,0},{0,7500,0}},out=s;const tinycowell_model m={398600441500000.,fail,NULL};CHECK(tinycowell_fixed(&m,0,1,&s,1,10,&out,w,54)==TINYCOWELL_FORCE_ERROR);CHECK(!memcmp(&out,&s,sizeof s));CHECK(tinycowell_fixed(&m,0,100,&s,1,1,&out,w,54)==TINYCOWELL_INTEGRATION_ERROR);CHECK(tinycowell_adaptive(&m,0,1,&s,&o,&out,&stats,w,1)==TINYCOWELL_INVALID);}
    printf("tinycowell: SciPy DOP853 adaptive maxima %.12g m %.12g m/s\n",maxp,maxv);return 0;}
