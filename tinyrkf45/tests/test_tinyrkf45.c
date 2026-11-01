/* SPDX-License-Identifier: MIT */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL %s:%d %s\n",__FILE__,__LINE__,#c); return 1; } } while (0)

#include "tinyrkf45/tinyrkf45.h"
#include "data/reference.h"
static int rhs(double t,const double *y,size_t n,void *ctx,double *dy){int model=*(int*)ctx;(void)n;(void)t;
    if(model==0){dy[0]=-.7*y[0];return 0;}
    if(model==1){dy[0]=y[1];dy[1]=-2.3*y[0];return 0;}
    if(model==2){dy[0]=10*(y[1]-y[0]);dy[1]=y[0]*(28-y[2])-y[1];dy[2]=y[0]*y[1]-(8.0/3)*y[2];return 0;}
    return 1;}
int main(void){double w[288],out[32],maximum=0;tinyrkf45_stats stats;tinyrkf45_options o={1e-13,1e-13,.01,1e-12,.1,100000};
    for(size_t i=0;i<sizeof cases/sizeof cases[0];++i){const reference_case *c=&cases[i];int model=c->model;
        CHECK(tinyrkf45_integrate(rhs,&model,c->n,0,c->time,c->initial,&o,out,&stats,w,288)==0);CHECK(stats.final_time==c->time);
        for(size_t j=0;j<c->n;++j){double e=fabs(out[j]-c->expected[j])/fmax(1,fabs(c->expected[j]));maximum=fmax(e,maximum);if(e>=2e-9)fprintf(stderr,"case %zu component %zu time %.17g expected %.17g got %.17g error %.12g\n",i,j,c->time,c->expected[j],out[j],e);CHECK(e<2e-9);}}
    {double y[2]={1,2},out2[2]={123,456};int model=1;tinyrkf45_stats sentinel={1,2,3,4,5},s=sentinel;tinyrkf45_options bad=o;bad.max_attempts=1;
    CHECK(tinyrkf45_integrate(rhs,&model,2,0,20,y,&bad,out2,&s,w,288)==TINYRKF45_STEP_LIMIT);CHECK(out2[0]==123&&out2[1]==456);CHECK(!memcmp(&s,&sentinel,sizeof s));
    model=3;CHECK(tinyrkf45_integrate(rhs,&model,2,0,.1,y,&o,out2,&s,w,288)==TINYRKF45_DERIVATIVE_ERROR);
    CHECK(tinyrkf45_integrate(rhs,&model,2,0,0,y,&o,out2,&s,w,288)==0&&out2[0]==1&&out2[1]==2&&s.evaluations==0);
    CHECK(tinyrkf45_integrate(rhs,&model,2,0,.1,y,&o,out2,&s,w,2)==TINYRKF45_INVALID);}
    printf("tinyrkf45: SciPy DOP853 max scaled error %.12g\n",maximum);return 0;}
