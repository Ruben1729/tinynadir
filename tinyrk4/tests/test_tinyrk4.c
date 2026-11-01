/* SPDX-License-Identifier: MIT */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL %s:%d %s\n",__FILE__,__LINE__,#c); return 1; } } while (0)

#include "tinyrk4/tinyrk4.h"
#include "data/reference.h"
static int rhs(double t,const double *y,size_t n,void *ctx,double *dy){
    int model=*(int*)ctx;(void)n;(void)t;
    if(model==0){dy[0]=-.7*y[0];return 0;}
    if(model==1){dy[0]=y[1];dy[1]=-2.3*y[0];return 0;}
    if(model==2){dy[0]=10*(y[1]-y[0]);dy[1]=y[0]*(28-y[2])-y[1];dy[2]=y[0]*y[1]-(8.0/3)*y[2];return 0;}
    return 1;
}
int main(void){
    double w[192],out[32],maxerror=0;
    for(size_t i=0;i<sizeof cases/sizeof cases[0];++i){const reference_case *c=&cases[i];int model=c->model;
        CHECK(tinyrk4_step(rhs,&model,c->n,0,c->time,c->initial,out,w,192)==0);
        for(size_t j=0;j<c->n;++j){double e=fabs(out[j]-c->expected[j])/fmax(1,fabs(c->expected[j]));maxerror=fmax(e,maxerror);CHECK(e<c->limit);}
    }
    {double y[2]={1,2},copy[2]={123,456};int model=3;
    CHECK(tinyrk4_step(rhs,&model,2,0,.1,y,copy,w,192)==TINYRK4_DERIVATIVE_ERROR&&copy[0]==123&&copy[1]==456);
    CHECK(tinyrk4_step(rhs,&model,33,0,.1,y,copy,w,192)==TINYRK4_INVALID);
    CHECK(tinyrk4_step(rhs,&model,2,0,.1,y,copy,w,1)==TINYRK4_INVALID);
    CHECK(tinyrk4_step(rhs,&model,2,0,NAN,y,copy,w,192)==TINYRK4_INVALID);
    CHECK(tinyrk4_step(rhs,&model,2,0,0,y,y,w,192)==0&&y[0]==1&&y[1]==2);}
    printf("tinyrk4: SciPy DOP853 comparison max scaled local error %.12g\n",maxerror);return 0;
}
