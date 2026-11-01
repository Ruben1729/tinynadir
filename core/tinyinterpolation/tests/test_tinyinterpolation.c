/* SPDX-License-Identifier: MIT */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL %s:%d %s\n",__FILE__,__LINE__,#c); return 1; } } while (0)
static int close_to(double a,double b,double tolerance) { return isfinite(a)&&fabs(a-b)<=tolerance*fmax(1.0,fabs(b)); }

#include "tinyinterpolation/tinyinterpolation.h"
#include "data/reference.h"
int main(void) {
    double value=123,deriv=456,x[3]={0,1,2},y[3]={0,1,4};
    for(size_t i=0;i<sizeof cases/sizeof cases[0];++i) {
        const reference_case *c=&cases[i];
        CHECK(tinyinterpolation_linear(c->x,c->y,4,c->query,&value)==TINYINTERPOLATION_OK);
        CHECK(close_to(value,c->linear,2e-12));
        CHECK(tinyinterpolation_lagrange(c->x,c->y,4,c->query,&value)==TINYINTERPOLATION_OK);
        CHECK(close_to(value,c->polynomial,2e-11));
        CHECK(tinyinterpolation_hermite(c->x[0],c->x[3],c->y[0],c->y[3],c->m0,c->m1,c->query,&value,&deriv)==TINYINTERPOLATION_OK);
        CHECK(close_to(value,c->hermite,2e-12));CHECK(close_to(deriv,c->derivative,2e-12));
    }
    CHECK(tinyinterpolation_linear(x,y,3,2,&value)==0&&value==4);
    value=123;CHECK(tinyinterpolation_linear(x,y,3,-1,&value)==TINYINTERPOLATION_RANGE&&value==123);
    x[1]=x[0];CHECK(tinyinterpolation_linear(x,y,3,0,&value)==TINYINTERPOLATION_INVALID&&value==123);
    x[1]=1;y[1]=NAN;CHECK(tinyinterpolation_lagrange(x,y,3,0,&value)==TINYINTERPOLATION_INVALID);
    CHECK(tinyinterpolation_hermite(0,0,0,1,0,0,0,&value,&deriv)==TINYINTERPOLATION_INVALID);
    CHECK(tinyinterpolation_hermite(0,1,0,1,0,0,NAN,&value,&deriv)==TINYINTERPOLATION_INVALID);
    CHECK(tinyinterpolation_linear(NULL,y,3,0,&value)==TINYINTERPOLATION_INVALID);
    CHECK(tinyinterpolation_linear(x,y,3,0,NULL)==TINYINTERPOLATION_INVALID);
    puts("tinyinterpolation: independent SciPy/NumPy fixtures and bounds pass");return 0;
}
