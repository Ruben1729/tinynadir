/* SPDX-License-Identifier: MIT */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL %s:%d %s\n",__FILE__,__LINE__,#c); return 1; } } while (0)

#include "tinygravity/tinygravity.h"
#include "data/reference.h"
#include "data/normalization.h"
static int normalization_tests(void){
    const unsigned degrees[]={0,1,2,8,16,24,32};
    double nc[561],ns[561],c[561],s[561],before_c[561],before_s[561];
    double maximum_relative=0;
    for(size_t i=0;i<561;++i){nc[i]=normalization_cases[i].c;ns[i]=normalization_cases[i].s;}
    for(size_t d=0;d<sizeof degrees/sizeof degrees[0];++d){
        size_t count=tinygravity_coefficient_count(degrees[d]);
        for(size_t i=0;i<561;++i){c[i]=123;s[i]=456;}
        CHECK(tinygravity_unnormalize(degrees[d],nc,ns,561,c,s)==TINYGRAVITY_OK);
        for(size_t i=0;i<count;++i){
            const normalization_case *r=&normalization_cases[i];
            double ec=fabs((c[i]-r->expected_c)/r->expected_c),es=fabs((s[i]-r->expected_s)/r->expected_s);
            CHECK(ec<2e-14&&es<2e-14);maximum_relative=fmax(maximum_relative,fmax(ec,es));
            /* Independently evaluated SciPy spherical harmonics/Legendre
             * bases check the angular field, not just the factorial formula. */
            CHECK(fabs((c[i]*r->unnormalized_basis-nc[i]*r->normalized_basis)/(nc[i]*r->normalized_basis))<3e-13);
            CHECK(fabs((s[i]*r->unnormalized_basis-ns[i]*r->normalized_basis)/(ns[i]*r->normalized_basis))<3e-13);
        }
        for(size_t i=count;i<561;++i)CHECK(c[i]==123&&s[i]==456);
    }
    memcpy(c,nc,sizeof c);memcpy(s,ns,sizeof s);
    CHECK(tinygravity_unnormalize(32,c,s,561,c,s)==TINYGRAVITY_OK);
    for(size_t i=0;i<561;++i){CHECK(fabs((c[i]-normalization_cases[i].expected_c)/normalization_cases[i].expected_c)<2e-14);CHECK(fabs((s[i]-normalization_cases[i].expected_s)/normalization_cases[i].expected_s)<2e-14);}
    memcpy(before_c,c,sizeof c);memcpy(before_s,s,sizeof s);
    /* Overflow late in the triangular table must leave both outputs intact. */
    nc[528]=DBL_MAX;
    CHECK(tinygravity_unnormalize(32,nc,ns,561,c,s)==TINYGRAVITY_ARITHMETIC);
    CHECK(!memcmp(c,before_c,sizeof c)&&!memcmp(s,before_s,sizeof s));
    nc[528]=normalization_cases[528].c;ns[528]=-DBL_MAX;
    CHECK(tinygravity_unnormalize(32,nc,ns,561,c,s)==TINYGRAVITY_ARITHMETIC);
    CHECK(!memcmp(c,before_c,sizeof c)&&!memcmp(s,before_s,sizeof s));
    memcpy(before_c,nc,sizeof nc);memcpy(before_s,ns,sizeof ns);
    CHECK(tinygravity_unnormalize(32,nc,ns,561,nc,ns)==TINYGRAVITY_ARITHMETIC);
    CHECK(!memcmp(nc,before_c,sizeof nc)&&!memcmp(ns,before_s,sizeof ns));
    ns[528]=normalization_cases[528].s;ns[560]=NAN;
    memcpy(before_c,c,sizeof c);memcpy(before_s,s,sizeof s);
    CHECK(tinygravity_unnormalize(32,nc,ns,561,c,s)==TINYGRAVITY_INVALID);
    CHECK(!memcmp(c,before_c,sizeof c)&&!memcmp(s,before_s,sizeof s));
    CHECK(tinygravity_unnormalize(33,nc,ns,561,c,s)==TINYGRAVITY_INVALID);
    CHECK(tinygravity_unnormalize(32,nc,ns,560,c,s)==TINYGRAVITY_INVALID);
    printf("tinygravity: 561 independent normalization / SciPy angular-basis cases through degree32; max relative %.3g\n",maximum_relative);
    return 0;
}
int main(void){double work[8*561],out[3],potential,maximum=0;
    CHECK(normalization_tests()==0);
    for(size_t i=0;i<sizeof cases/sizeof cases[0];++i){const reference_case *r=&cases[i];const reference_model *m=&models[r->model];CHECK(tinygravity_evaluate(398600441500000.0,6378136.3,m->degree,m->degree,m->c,m->s,m->count,r->r,out,&potential,work,8*561)==0);for(size_t j=0;j<3;++j){double e=fabs(out[j]-r->acceleration[j]);maximum=fmax(maximum,e);if(e>=3e-11)fprintf(stderr,"case %zu degree %u j%zu got %.17g expected %.17g error %.12g\n",i,m->degree,j,out[j],r->acceleration[j],e);CHECK(e<3e-11);}CHECK(potential>0);}
    {double c[6]={1,0,0,-1.08262668e-3,0,0},s[6]={0},r[3]={0,0,7000000};double mu=398600441500000.,a=6378136.3,radius=r[2],expected=-mu/(radius*radius)*(1-3*1.08262668e-3*(a/radius)*(a/radius));CHECK(tinygravity_evaluate(mu,a,2,0,c,s,6,r,out,&potential,work,8*561)==0);CHECK(fabs(out[2]-expected)<1e-13&&out[0]==0&&out[1]==0);out[0]=123;r[2]=0;CHECK(tinygravity_evaluate(mu,a,2,0,c,s,6,r,out,NULL,work,8*561)==TINYGRAVITY_SINGULAR&&out[0]==123);CHECK(tinygravity_evaluate(mu,a,2,0,c,s,6,r,out,NULL,work,1)==TINYGRAVITY_INVALID);}
    {double c[6]={1,0,0,-1.08262668e-3/sqrt(5),0,0},s[6]={0};CHECK(tinygravity_unnormalize(2,c,s,6,c,s)==0);CHECK(fabs(c[3]+1.08262668e-3)<1e-18);}
    printf("tinygravity: independent Vallado Montenbruck comparison max acceleration error %.12g m/s2\n",maximum);return 0;}
