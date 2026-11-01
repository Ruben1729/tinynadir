/* SPDX-License-Identifier: MIT */
#include "tinybatchls/tinybatchls.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "data/reference.h"
#define CHECK(x) do { if (!(x)) {fprintf(stderr,"line %d: %s\n",__LINE__,#x);return 1;} } while(0)
typedef struct {unsigned kind;const double *x;} context;
static int model(const double *p,size_t n,double *y,size_t m,void *data)
{context *c=data;size_t i;(void)n;for(i=0;i<m;++i){double x=c->x[i];if(c->kind==0)y[i]=p[0]*x+p[1];else if(c->kind==1)y[i]=p[0]*exp(-p[1]*x)+p[2];else if(c->kind==2)y[i]=p[0]/(1+exp(-p[1]*(x-p[2])));else if(c->kind==3)y[i]=p[0]+p[1];else return 0;}return 1;}
int main(void)
{
    size_t i,j;double maximum=0,max_cov=0;tinybatchls_workspace w;
    for(i=0;i<sizeof reference/sizeof reference[0];++i){tinybatchls_options o;tinybatchls_result out;context c={reference[i].kind,reference[i].x};tinybatchls_status status;CHECK(tinybatchls_default_options(reference[i].n,reference[i].m,&o)==TINYBATCHLS_OK);o.step_tolerance=1e-11;o.gradient_tolerance=1e-6;status=tinybatchls_solve(&o,model,&c,reference[i].y,reference[i].sigma,reference[i].initial,&w,&out);if(status!=TINYBATCHLS_OK)fprintf(stderr,"fit %zu status %d\n",i,(int)status);CHECK(status==TINYBATCHLS_OK);for(j=0;j<reference[i].n;++j){double e=fabs(out.parameters[j]-reference[i].parameters[j]);if(e>maximum)maximum=e;CHECK(e<2e-7);}for(j=0;j<reference[i].n*reference[i].n;++j){double scale=sqrt(fabs(reference[i].covariance[(j/reference[i].n)*reference[i].n+j/reference[i].n]*reference[i].covariance[(j%reference[i].n)*reference[i].n+j%reference[i].n]));double e=fabs(out.covariance[j]-reference[i].covariance[j])/scale;if(e>max_cov)max_cov=e;CHECK(e<3e-6);}CHECK(fabs(out.weighted_sum_squares-reference[i].cost)<2e-8);CHECK(out.iterations<=o.max_iterations);CHECK(out.model_evaluations<=1+(2*o.parameter_count+1)*o.max_iterations+2*o.parameter_count);}
    {tinybatchls_options o;tinybatchls_result out={0},saved=out;context c={3,reference[0].x};CHECK(tinybatchls_default_options(2,9,&o)==TINYBATCHLS_OK);CHECK(tinybatchls_solve(&o,model,&c,reference[0].y,reference[0].sigma,reference[0].initial,&w,&out)==TINYBATCHLS_RANK_DEFICIENT);CHECK(!memcmp(&out,&saved,sizeof out));c.kind=4;CHECK(tinybatchls_solve(&o,model,&c,reference[0].y,reference[0].sigma,reference[0].initial,&w,&out)==TINYBATCHLS_MODEL_ERROR);CHECK(!memcmp(&out,&saved,sizeof out));}
    printf("Batch LS: SciPy linear/exponential/logistic fits; max parameter %.9g, normalized covariance %.9g; workspace %zu bytes\n",maximum,max_cov,sizeof w);return 0;
}
