/* SPDX-License-Identifier: MIT */
#include "tinybatchls/tinybatchls.h"
#include <math.h>
#include <float.h>
static int solve(size_t n,double *A,double *b)
{
    double scale=0;size_t i,j,k;
    for(i=0;i<n*n;++i)if(fabs(A[i])>scale)scale=fabs(A[i]);
    if(!isfinite(scale)||scale==0)return 0;
    for(k=0;k<n;++k){size_t pivot=k;double best=fabs(A[k*n+k]);for(i=k+1;i<n;++i)if(fabs(A[i*n+k])>best){best=fabs(A[i*n+k]);pivot=i;}if(best<=1e-14*scale)return 0;if(pivot!=k){double temp=b[k];b[k]=b[pivot];b[pivot]=temp;for(j=k;j<n;++j){temp=A[k*n+j];A[k*n+j]=A[pivot*n+j];A[pivot*n+j]=temp;}}for(i=k+1;i<n;++i){double f=A[i*n+k]/A[k*n+k];A[i*n+k]=0;for(j=k+1;j<n;++j)A[i*n+j]-=f*A[k*n+j];b[i]-=f*b[k];}}
    for(i=n;i-- >0;){double value=b[i];for(j=i+1;j<n;++j)value-=A[i*n+j]*b[j];b[i]=value/A[i*n+i];if(!isfinite(b[i]))return 0;}
    return 1;
}
tinybatchls_status tinybatchls_default_options(size_t n,size_t m,tinybatchls_options *out)
{
    tinybatchls_options o={0};size_t i;if(!out||n<1||n>8||m<n||m>64)return TINYBATCHLS_INVALID;
    o.parameter_count=n;o.observation_count=m;o.max_iterations=100;
    for(i=0;i<8;++i)o.parameter_scale[i]=1;
    o.relative_difference_step=1e-5;o.step_tolerance=1e-10;o.gradient_tolerance=1e-7;o.cost_tolerance=1e-12;o.initial_damping=1e-3;
    *out=o;return TINYBATCHLS_OK;
}
static int predict(const tinybatchls_options *o,tinybatchls_model model,void *context,const double *x,double *y,size_t *count)
{
    size_t i;for(i=0;i<o->observation_count;++i)y[i]=NAN;++*count;
    if(!model(x,o->parameter_count,y,o->observation_count,context))return 0;
    for(i=0;i<o->observation_count;++i)if(!isfinite(y[i]))return 0;
    return 1;
}
static int cost(const tinybatchls_options *o,const double *observations,const double *sigma,const double *prediction,double *sum)
{
    size_t i;double c=0;for(i=0;i<o->observation_count;++i){double r=(observations[i]-prediction[i])/sigma[i];c+=r*r;}if(!isfinite(c))return 0;*sum=c;return 1;
}
static tinybatchls_status linearize(const tinybatchls_options *o,tinybatchls_model model,void *context,const double *x,const double *y,const double *sigma,tinybatchls_workspace *w,size_t *evaluations)
{
    size_t n=o->parameter_count,m=o->observation_count,i,j,k;double trial[8];
    for(j=0;j<n;++j){double h=o->relative_difference_step*fmax(fabs(x[j]),o->parameter_scale[j]);if(!isfinite(h)||!(h>0)||!isfinite(x[j]+h)||!isfinite(x[j]-h)||x[j]+h==x[j]||x[j]-h==x[j])return TINYBATCHLS_INVALID;for(i=0;i<n;++i)trial[i]=x[i];trial[j]=x[j]+h;if(!predict(o,model,context,trial,w->plus,evaluations))return TINYBATCHLS_MODEL_ERROR;trial[j]=x[j]-h;if(!predict(o,model,context,trial,w->minus,evaluations))return TINYBATCHLS_MODEL_ERROR;for(i=0;i<m;++i){w->jacobian[n*i+j]=((w->plus[i]-w->minus[i])/(2*h))*o->parameter_scale[j]/sigma[i];if(!isfinite(w->jacobian[n*i+j]))return TINYBATCHLS_INVALID;}}
    for(j=0;j<n;++j){w->gradient[j]=0;for(k=0;k<n;++k)w->normal[n*j+k]=0;}
    for(i=0;i<m;++i){double residual=(y[i]-w->prediction[i])/sigma[i];for(j=0;j<n;++j){double J=w->jacobian[n*i+j];w->gradient[j]+=J*residual;for(k=0;k<n;++k)w->normal[n*j+k]+=J*w->jacobian[n*i+k];}}
    for(j=0;j<n;++j){if(!isfinite(w->gradient[j]))return TINYBATCHLS_INVALID;for(k=0;k<n;++k)if(!isfinite(w->normal[n*j+k]))return TINYBATCHLS_INVALID;}
    return TINYBATCHLS_OK;
}
tinybatchls_status tinybatchls_solve(const tinybatchls_options *o,tinybatchls_model model,void *context,const double *y,const double *sigma,const double *initial,tinybatchls_workspace *w,tinybatchls_result *out)
{
    tinybatchls_result result={0};size_t n,m,i,j,iteration,evaluations=0;double x[8],chi,lambda;int converged=0;tinybatchls_status status;
    if(!o||!model||!y||!sigma||!initial||!w||!out)return TINYBATCHLS_INVALID;
    n=o->parameter_count;m=o->observation_count;
    if(n<1||n>8||m<n||m>64||o->max_iterations<1||o->max_iterations>100||
       !isfinite(o->relative_difference_step)||o->relative_difference_step<1e-10||o->relative_difference_step>.01||
       !isfinite(o->step_tolerance)||o->step_tolerance<=0||!isfinite(o->gradient_tolerance)||o->gradient_tolerance<=0||
       !isfinite(o->cost_tolerance)||o->cost_tolerance<=0||!isfinite(o->initial_damping)||o->initial_damping<0||o->initial_damping>1e12)return TINYBATCHLS_INVALID;
    for(i=0;i<n;++i){if(!isfinite(initial[i])||!isfinite(o->parameter_scale[i])||o->parameter_scale[i]<=0||!isfinite(initial[i]/o->parameter_scale[i]))return TINYBATCHLS_INVALID;x[i]=initial[i];}
    for(i=0;i<m;++i)if(!isfinite(y[i])||!isfinite(sigma[i])||sigma[i]<=0)return TINYBATCHLS_INVALID;
    if(!predict(o,model,context,x,w->prediction,&evaluations))return TINYBATCHLS_MODEL_ERROR;
    if(!cost(o,y,sigma,w->prediction,&chi))return TINYBATCHLS_INVALID;
    lambda=o->initial_damping;
    for(iteration=0;iteration<o->max_iterations;++iteration){double gradient=0,step=0,xnorm=0,trial_cost;
        status=linearize(o,model,context,x,y,sigma,w,&evaluations);if(status!=TINYBATCHLS_OK)return status;
        for(i=0;i<n;++i){gradient=fmax(gradient,fabs(w->gradient[i]));xnorm=fmax(xnorm,fabs(x[i]/o->parameter_scale[i]));}
        if(gradient<=o->gradient_tolerance){converged=1;break;}
        for(i=0;i<n*n;++i)w->factor[i]=w->normal[i];
        for(i=0;i<n;++i){if(!(w->normal[n*i+i]>0))return TINYBATCHLS_RANK_DEFICIENT;w->factor[n*i+i]+=lambda*w->normal[n*i+i];w->delta[i]=w->gradient[i];}
        if(!solve(n,w->factor,w->delta))return TINYBATCHLS_RANK_DEFICIENT;
        for(i=0;i<n;++i){step=fmax(step,fabs(w->delta[i]));w->trial[i]=x[i]+o->parameter_scale[i]*w->delta[i];if(!isfinite(w->trial[i])||!isfinite(w->trial[i]/o->parameter_scale[i]))return TINYBATCHLS_INVALID;}
        if(step<=o->step_tolerance*(1+xnorm) && lambda<=1e-3){converged=1;break;}
        if(!predict(o,model,context,w->trial,w->trial_prediction,&evaluations)) {lambda=lambda==0 ? 1e-3 : lambda*10;if(lambda>1e12)return TINYBATCHLS_NO_CONVERGENCE;continue;}
        if(!cost(o,y,sigma,w->trial_prediction,&trial_cost))return TINYBATCHLS_INVALID;
        if(trial_cost<chi){double improvement=chi-trial_cost;for(i=0;i<n;++i)x[i]=w->trial[i];for(i=0;i<m;++i)w->prediction[i]=w->trial_prediction[i];chi=trial_cost;lambda=fmax(1e-15,lambda*.3);if(improvement<=o->cost_tolerance*(1+chi)&&step<=sqrt(o->step_tolerance)*(1+xnorm)){converged=1;break;}}
        else {lambda=lambda==0 ? 1e-3 : lambda*10;if(lambda>1e12)return TINYBATCHLS_NO_CONVERGENCE;}
    }
    if(!converged)return TINYBATCHLS_NO_CONVERGENCE;
    status=linearize(o,model,context,x,y,sigma,w,&evaluations);if(status!=TINYBATCHLS_OK)return status;
    for(j=0;j<n;++j){for(i=0;i<n*n;++i)w->factor[i]=w->normal[i];for(i=0;i<n;++i)w->delta[i]=i==j ? 1 : 0;if(!solve(n,w->factor,w->delta))return TINYBATCHLS_RANK_DEFICIENT;for(i=0;i<n;++i){result.covariance[n*i+j]=(o->parameter_scale[i]*w->delta[i])*o->parameter_scale[j];if(!isfinite(result.covariance[n*i+j]))return TINYBATCHLS_INVALID;}}
    for(i=0;i<n;++i)result.parameters[i]=x[i];
    result.weighted_sum_squares=chi;result.iterations=iteration+1;result.model_evaluations=evaluations;
    *out=result;return TINYBATCHLS_OK;
}
