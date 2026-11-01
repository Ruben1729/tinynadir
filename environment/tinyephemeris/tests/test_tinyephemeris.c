/* SPDX-License-Identifier: MIT */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL %s:%d %s\n",__FILE__,__LINE__,#c); return 1; } } while (0)
static int close_to(double a,double b,double tolerance) { return isfinite(a)&&fabs(a-b)<=tolerance*fmax(1.0,fabs(b)); }

#include "tinyephemeris/tinyephemeris.h"
#include "data/reference.h"
int main(void){double p[3],v[3],maximum=0,maxvelocity=0;
    for(size_t i=0;i<sizeof cases/sizeof cases[0];++i){const reference_case *c=&cases[i];CHECK(tinyephemeris_evaluate(c->samples,4,c->time,p,v)==0);for(size_t j=0;j<3;++j){CHECK(close_to(p[j],c->p[j],3e-13));CHECK(fabs(v[j]-c->v[j])<1e-9+3e-13*fabs(c->v[j]));maximum=fmax(maximum,fabs(p[j]-c->p[j]));maxvelocity=fmax(maxvelocity,fabs(v[j]-c->v[j]));}}
    {tinyephemeris_sample samples[2]={{0,{0,0,0},{1,2,3}},{1,{1,2,3},{1,2,3}}};CHECK(tinyephemeris_evaluate(samples,2,0,p,v)==0&&p[0]==0&&v[2]==3);p[0]=123;CHECK(tinyephemeris_evaluate(samples,2,-1,p,v)==TINYEPHEMERIS_RANGE&&p[0]==123);samples[1].time_s=0;CHECK(tinyephemeris_evaluate(samples,2,0,p,v)==TINYEPHEMERIS_INVALID);}
    printf("tinyephemeris: SciPy Hermite ephemeris max position/velocity differences %.12g m %.12g m/s\n",maximum,maxvelocity);return 0;}
