/* SPDX-License-Identifier: MIT */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL %s:%d %s\n",__FILE__,__LINE__,#c); return 1; } } while (0)
static int close_to(double a,double b,double tolerance) { return isfinite(a)&&fabs(a-b)<=tolerance*fmax(1.0,fabs(b)); }

#include "tinyspaceweather/tinyspaceweather.h"
#include "data/reference.h"
int main(void){tinyspaceweather_sample out;double value;
    for(size_t i=0;i<sizeof cases/sizeof cases[0];++i){const reference_case *c=&cases[i];CHECK(tinyspaceweather_interpolate(c->samples,4,c->time,&out)==0);CHECK(close_to(out.f107_sfu,c->expected[0],3e-14));CHECK(close_to(out.f107_average_sfu,c->expected[1],3e-14));CHECK(close_to(out.ap,c->expected[2],3e-14));}
    for(size_t i=0;i<sizeof conversions/sizeof conversions[0];++i){CHECK(tinyspaceweather_kp_to_ap(conversions[i][0],&value)==0);CHECK(close_to(value,conversions[i][1],3e-14));CHECK(tinyspaceweather_ap_to_kp(conversions[i][1],&value)==0);CHECK(close_to(value,conversions[i][0],3e-14));}
    value=123;CHECK(tinyspaceweather_kp_to_ap(10,&value)==TINYSPACEWEATHER_RANGE&&value==123);CHECK(tinyspaceweather_ap_to_kp(-1,&value)==TINYSPACEWEATHER_RANGE);CHECK(tinyspaceweather_kp_to_ap(NAN,&value)==TINYSPACEWEATHER_INVALID);
    out.time_s=123;CHECK(tinyspaceweather_interpolate(cases[0].samples,4,-10000,&out)==TINYSPACEWEATHER_RANGE&&out.time_s==123);
    puts("tinyspaceweather: NumPy table references and Kp/Ap conversions pass");return 0;}
