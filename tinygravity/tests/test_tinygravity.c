/* SPDX-License-Identifier: MIT */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL %s:%d %s\n",__FILE__,__LINE__,#c); return 1; } } while (0)

#include "tinygravity/tinygravity.h"
#include "data/reference.h"
int main(void){double work[8*561],out[3],potential,maximum=0;
    for(size_t i=0;i<sizeof cases/sizeof cases[0];++i){const reference_case *r=&cases[i];const reference_model *m=&models[r->model];CHECK(tinygravity_evaluate(398600441500000.0,6378136.3,m->degree,m->degree,m->c,m->s,m->count,r->r,out,&potential,work,8*561)==0);for(size_t j=0;j<3;++j){double e=fabs(out[j]-r->acceleration[j]);maximum=fmax(maximum,e);if(e>=3e-11)fprintf(stderr,"case %zu degree %u j%zu got %.17g expected %.17g error %.12g\n",i,m->degree,j,out[j],r->acceleration[j],e);CHECK(e<3e-11);}CHECK(potential>0);}
    {double c[6]={1,0,0,-1.08262668e-3,0,0},s[6]={0},r[3]={0,0,7000000};double mu=398600441500000.,a=6378136.3,radius=r[2],expected=-mu/(radius*radius)*(1-3*1.08262668e-3*(a/radius)*(a/radius));CHECK(tinygravity_evaluate(mu,a,2,0,c,s,6,r,out,&potential,work,8*561)==0);CHECK(fabs(out[2]-expected)<1e-13&&out[0]==0&&out[1]==0);out[0]=123;r[2]=0;CHECK(tinygravity_evaluate(mu,a,2,0,c,s,6,r,out,NULL,work,8*561)==TINYGRAVITY_SINGULAR&&out[0]==123);CHECK(tinygravity_evaluate(mu,a,2,0,c,s,6,r,out,NULL,work,1)==TINYGRAVITY_INVALID);}
    {double c[6]={1,0,0,-1.08262668e-3/sqrt(5),0,0},s[6]={0};CHECK(tinygravity_unnormalize(2,c,s,6,c,s)==0);CHECK(fabs(c[3]+1.08262668e-3)<1e-18);}
    printf("tinygravity: independent Vallado Montenbruck comparison max acceleration error %.12g m/s2\n",maximum);return 0;}
