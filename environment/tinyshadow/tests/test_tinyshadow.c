/* SPDX-License-Identifier: MIT */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL %s:%d %s\n",__FILE__,__LINE__,#c); return 1; } } while (0)

#include "tinyshadow/tinyshadow.h"
#include "data/reference.h"
static int light(double t,void *ctx,double *value){(void)ctx;*value=.5+.5*sin(t);return 0;}
int main(void){double out,maximum=0;int blocked;
    for(size_t i=0;i<sizeof cases/sizeof cases[0];++i){const reference_case *c=&cases[i];CHECK(tinyshadow_sun_fraction(c->r,c->sun,c->rb,c->rs,&out)==0);maximum=fmax(maximum,fabs(out-c->fraction));CHECK(fabs(out-c->fraction)<5e-8);}
    for(size_t i=0;i<sizeof los/sizeof los[0];++i){CHECK(tinyshadow_segment_ellipsoid(los[i].a,los[i].b,los[i].axes,&blocked)==0);CHECK(blocked==los[i].blocked);}
    CHECK(tinyshadow_event(light,NULL,-1,1,.7,1e-10,64,&out)==0);CHECK(fabs(out-event_reference)<1e-10);
    out=123;CHECK(tinyshadow_event(light,NULL,0,.1,.9,1e-10,64,&out)==TINYSHADOW_UNBRACKETED&&out==123);CHECK(tinyshadow_event(light,NULL,-1,1,.7,1e-15,1,&out)==TINYSHADOW_NO_CONVERGENCE);
    {double r[3]={-7000000,0,0},sun[3]={1,0,0};CHECK(tinyshadow_cylindrical(r,sun,6378136,&blocked)==0&&blocked==1);r[0]=-r[0];CHECK(tinyshadow_cylindrical(r,sun,6378136,&blocked)==0&&blocked==0);out=123;r[0]=0;CHECK(tinyshadow_sun_fraction(r,sun,6378136,695700000,&out)==TINYSHADOW_GEOMETRY&&out==123);}
    printf("tinyshadow: independent disc integration maximum fraction error %.12g\n",maximum);return 0;}
