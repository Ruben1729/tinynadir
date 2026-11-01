/* SPDX-License-Identifier: MIT */
#include "tinysgp4fit/tinysgp4fit.h"
#include <math.h>
#include <float.h>
typedef struct {
    tinytle_elements seed;
    const tinysgp4fit_observation *observations;
    size_t count;
    int velocity,bstar;
    tinysgp4fit_workspace *workspace;
} fit_context;
static double degrees(double x){x=fmod(x,360);return x<0 ? x+360 : x;}
static int decode(const double *x,const fit_context *context,tinytle_elements *out)
{
    tinytle_elements e=context->seed;double temp;
    if(x[0]<=1e-6||x[0]>=180-1e-6||x[5]>log(DBL_MAX))return 0;
    if(x[2]>=0)e.eccentricity=1/(1+exp(-x[2]));else {temp=exp(x[2]);e.eccentricity=temp/(1+temp);}
    e.inclination_deg=x[0];e.raan_deg=degrees(x[1]);e.arg_perigee_deg=degrees(x[3]);e.mean_anomaly_deg=degrees(x[4]);
    e.mean_motion_rev_day=exp(x[5]);if(context->bstar)e.bstar=x[6];
    if(e.eccentricity<=1e-12||e.eccentricity>=1-1e-12||!isfinite(e.mean_motion_rev_day)||!(e.mean_motion_rev_day>0))return 0;
    *out=e;return 1;
}
static int model(const double *x,size_t n,double *prediction,size_t m,void *data)
{
    fit_context *context=data;tinytle_elements e;size_t i,j,k=0;size_t expected=context->velocity ? 6*context->count : 3*context->count;
    if(n!=(context->bstar?7u:6u)||m!=expected||!decode(x,context,&e)||tinytle_initialize_elements(&e,&context->workspace->satellite)!=TINYTLE_OK)return 0;
    for(i=0;i<context->count;++i){if(tgp4_propagate(&context->workspace->satellite,context->observations[i].minutes_since_epoch,&context->workspace->state)!=TGP4_OK)return 0;for(j=0;j<3;++j)prediction[k++]=context->workspace->state.position_km[j]*1000;if(context->velocity)for(j=0;j<3;++j)prediction[k++]=context->workspace->state.velocity_km_s[j]*1000;}
    return 1;
}
tinysgp4fit_status tinysgp4fit_default_options(tinysgp4fit_options *out)
{
    tinysgp4fit_options o={100,0,1,1e-3};if(!out)return TINYSGP4FIT_INVALID;*out=o;return TINYSGP4FIT_OK;
}
tinysgp4fit_status tinysgp4fit_fit(const tinytle_elements *seed,const tinysgp4fit_observation *observations,size_t count,const tinysgp4fit_options *options,tinysgp4fit_workspace *w,tinysgp4fit_result *out)
{
    tinybatchls_options batch;tinybatchls_result fit;fit_context context;tinysgp4fit_result result={0};
    double initial[8]={0},derivative[7],tmin,tmax;size_t n,m,i,j,k=0;tinybatchls_status status;
    if(!seed||!observations||!options||!w||!out||count<2||count>10||options->max_iterations<1||options->max_iterations>100||
       (options->fit_bstar!=0&&options->fit_bstar!=1)||(options->use_velocity!=0&&options->use_velocity!=1)||
       (!options->use_velocity&&count<3)||!isfinite(options->initial_damping)||options->initial_damping<0||options->initial_damping>1e12)return TINYSGP4FIT_INVALID;
    if(!isfinite(seed->eccentricity)||seed->eccentricity<=1e-12||seed->eccentricity>=1-1e-12||!isfinite(seed->inclination_deg)||seed->inclination_deg<=1e-6||seed->inclination_deg>=180-1e-6)return TINYSGP4FIT_UNOBSERVABLE;
    if(tinytle_initialize_elements(seed,&w->satellite)!=TINYTLE_OK)return TINYSGP4FIT_INVALID;
    n=options->fit_bstar?7:6;m=count*(options->use_velocity?6:3);
    if(tinybatchls_default_options(n,m,&batch)!=TINYBATCHLS_OK)return TINYSGP4FIT_INVALID;
    batch.max_iterations=options->max_iterations;batch.initial_damping=options->initial_damping;
    batch.relative_difference_step=1e-6;batch.step_tolerance=1e-13;batch.gradient_tolerance=1e-5;batch.cost_tolerance=1e-8;
    batch.parameter_scale[0]=10;batch.parameter_scale[1]=10;batch.parameter_scale[2]=1;
    batch.parameter_scale[3]=10;batch.parameter_scale[4]=10;batch.parameter_scale[5]=.01;batch.parameter_scale[6]=.1;
    initial[0]=seed->inclination_deg;initial[1]=seed->raan_deg;initial[2]=log(seed->eccentricity/(1-seed->eccentricity));
    initial[3]=seed->arg_perigee_deg;initial[4]=seed->mean_anomaly_deg;initial[5]=log(seed->mean_motion_rev_day);initial[6]=seed->bstar;
    tmin=tmax=observations[0].minutes_since_epoch;
    for(i=0;i<count;++i){double t=observations[i].minutes_since_epoch;if(!isfinite(t)||fabs(t)>10080)return TINYSGP4FIT_INVALID;tmin=fmin(tmin,t);tmax=fmax(tmax,t);
        for(j=0;j<3;++j){if(!isfinite(observations[i].position_m[j])||!isfinite(observations[i].position_sigma_m[j])||observations[i].position_sigma_m[j]<=0)return TINYSGP4FIT_INVALID;w->scalar_observations[k]=observations[i].position_m[j];w->sigma[k++]=observations[i].position_sigma_m[j];}
        if(options->use_velocity)for(j=0;j<3;++j){if(!isfinite(observations[i].velocity_m_s[j])||!isfinite(observations[i].velocity_sigma_m_s[j])||observations[i].velocity_sigma_m_s[j]<=0)return TINYSGP4FIT_INVALID;w->scalar_observations[k]=observations[i].velocity_m_s[j];w->sigma[k++]=observations[i].velocity_sigma_m_s[j];}
    }
    if(options->fit_bstar&&tmax==tmin)return TINYSGP4FIT_UNOBSERVABLE;
    context.seed=*seed;context.observations=observations;context.count=count;context.velocity=options->use_velocity;context.bstar=options->fit_bstar;context.workspace=w;
    status=tinybatchls_solve(&batch,model,&context,w->scalar_observations,w->sigma,initial,&w->batch,&fit);
    if(status==TINYBATCHLS_RANK_DEFICIENT)return TINYSGP4FIT_UNOBSERVABLE;
    if(status==TINYBATCHLS_MODEL_ERROR)return TINYSGP4FIT_MODEL_ERROR;
    if(status==TINYBATCHLS_NO_CONVERGENCE)return TINYSGP4FIT_NO_CONVERGENCE;
    if(status!=TINYBATCHLS_OK)return TINYSGP4FIT_INVALID;
    if(!decode(fit.parameters,&context,&result.elements))return TINYSGP4FIT_MODEL_ERROR;
    for(i=0;i<7;++i)derivative[i]=1;
    derivative[2]=result.elements.eccentricity*(1-result.elements.eccentricity);derivative[5]=result.elements.mean_motion_rev_day;
    for(i=0;i<n;++i)for(j=0;j<n;++j){result.covariance[n*i+j]=(derivative[i]*fit.covariance[n*i+j])*derivative[j];if(!isfinite(result.covariance[n*i+j]))return TINYSGP4FIT_INVALID;}
    result.weighted_sum_squares=fit.weighted_sum_squares;result.normalized_rms=sqrt(fit.weighted_sum_squares/(double)m);
    result.parameter_count=n;result.scalar_observation_count=m;result.iterations=fit.iterations;result.model_evaluations=fit.model_evaluations;
    *out=result;return TINYSGP4FIT_OK;
}
