/* SPDX-License-Identifier: MIT */
#include "tinymekf/tinymekf.h"
#include "tinylinalg/tinylinalg.h"
#include "tinyrotation/tinyrotation.h"
#include "tinygyro/tinygyro.h"
#include <math.h>
#include <stddef.h>
#include <string.h>

static int finite_array(const double *a,size_t n){if(!a)return 0;for(size_t i=0;i<n;++i)if(!isfinite(a[i]))return 0;return 1;}
static void identity(double a[36]){for(size_t i=0;i<36;++i)a[i]=0;for(size_t i=0;i<6;++i)a[i*6+i]=1;}
static void multiply(const double a[36],const double b[36],double out[36],int transpose_b){
    for(size_t i=0;i<6;++i)for(size_t j=0;j<6;++j){double x=0;for(size_t k=0;k<6;++k)x+=a[i*6+k]*(transpose_b?b[j*6+k]:b[k*6+j]);out[i*6+j]=x;}
}
static void skew(const double v[3],double out[9]){
    out[0]=0;out[1]=-v[2];out[2]=v[1];out[3]=v[2];out[4]=0;out[5]=-v[0];out[6]=-v[1];out[7]=v[0];out[8]=0;
}
static int covariance_valid(const double p[36],tinymekf_workspace *work){
    return tinylinalg_cholesky(6,p,work->scratch,work->scratch+36,36)==TINYLINALG_OK;
}
static int symmetrize(double p[36]){
    if(!finite_array(p,36))return 0;
    for(size_t i=0;i<6;++i)for(size_t j=i+1;j<6;++j)p[i*6+j]=p[j*6+i]=.5*p[i*6+j]+.5*p[j*6+i];
    return 1;
}
static tinymekf_status check_state(const tinymekf_state *state,tinymekf_workspace *work,double normalized[4]){
    if(!state||!work||!finite_array(state->bias,3)||tinyrotation_normalize(state->q,normalized)!=TINYROTATION_OK)return TINYMEKF_INVALID_ARGUMENT;
    if(!covariance_valid(state->covariance,work))return TINYMEKF_INVALID_COVARIANCE;
    return TINYMEKF_OK;
}
tinymekf_status tinymekf_initialize(tinymekf_state *state,const double q[4],const double bias[3],const double covariance[36],tinymekf_workspace *work){
    tinymekf_state next;
    if(!state||!work||!finite_array(bias,3)||tinyrotation_normalize(q,next.q)!=TINYROTATION_OK)return TINYMEKF_INVALID_ARGUMENT;
    if(!covariance_valid(covariance,work))return TINYMEKF_INVALID_COVARIANCE;
    memcpy(next.bias,bias,sizeof next.bias);memcpy(next.covariance,covariance,sizeof next.covariance);*state=next;return TINYMEKF_OK;
}
tinymekf_status tinymekf_predict(tinymekf_state *state,const double gyro[3],double dt,
    const double gyro_psd[3],const double bias_psd[3],tinymekf_workspace *work){
    tinymekf_state next;double normalized[4],omega[3],cross[9];tinymekf_status status;
    if(!state||!work||!isfinite(dt)||dt<=0||dt>1||!finite_array(gyro,3)||!finite_array(gyro_psd,3)||!finite_array(bias_psd,3))return TINYMEKF_INVALID_ARGUMENT;
    for(size_t i=0;i<3;++i)if(gyro_psd[i]<0||bias_psd[i]<0)return TINYMEKF_INVALID_ARGUMENT;
    status=check_state(state,work,normalized);if(status!=TINYMEKF_OK)return status;
    next=*state;
    if(tinygyro_step(normalized,gyro,state->bias,dt,next.q)!=TINYGYRO_OK)return TINYMEKF_ARITHMETIC_ERROR;
    for(size_t i=0;i<3;++i){omega[i]=gyro[i]-state->bias[i];if(!isfinite(omega[i]))return TINYMEKF_ARITHMETIC_ERROR;}
    skew(omega,cross);identity(work->a);
    for(size_t i=0;i<3;++i){for(size_t j=0;j<3;++j)work->a[i*6+j]-=cross[i*3+j]*dt;work->a[i*6+i+3]=-dt;}
    multiply(work->a,state->covariance,work->b,0);multiply(work->b,work->a,work->c,1);
    for(size_t i=0;i<3;++i){work->c[i*6+i]+=gyro_psd[i]*dt;work->c[(i+3)*6+i+3]+=bias_psd[i]*dt;}
    if(!symmetrize(work->c)||!covariance_valid(work->c,work))return TINYMEKF_ARITHMETIC_ERROR;
    memcpy(next.covariance,work->c,sizeof next.covariance);*state=next;return TINYMEKF_OK;
}
tinymekf_status tinymekf_update_vector(tinymekf_state *state,const double reference[3],const double body[3],
    double variance,double gate,double *nis,tinymekf_workspace *work){
    tinymekf_state next;double normalized[4],r[3],b[3],predicted[3],cross[9],residual[3],solved[3],correction[6]={0},statistic=0,delta[4],rotation[3];tinymekf_status status;
    if(!state||!work||!finite_array(reference,3)||!finite_array(body,3)||!isfinite(variance)||variance<=0||!isfinite(gate)||gate<0)return TINYMEKF_INVALID_ARGUMENT;
    status=check_state(state,work,normalized);if(status!=TINYMEKF_OK)return status;
    if(tinylinalg_unit3(reference,r)!=TINYLINALG_OK||tinylinalg_unit3(body,b)!=TINYLINALG_OK)return TINYMEKF_DEGENERATE;
    if(tinyrotation_apply(normalized,r,predicted)!=TINYROTATION_OK)return TINYMEKF_ARITHMETIC_ERROR;
    skew(predicted,cross);
    for(size_t i=0;i<3;++i){residual[i]=b[i]-predicted[i];for(size_t j=0;j<6;++j)work->measurement[i*6+j]=j<3?cross[i*3+j]:0;}
    for(size_t i=0;i<3;++i)for(size_t j=0;j<3;++j){
        double value=i==j?variance:0;
        for(size_t k=0;k<6;++k)for(size_t l=0;l<6;++l)value+=work->measurement[i*6+k]*state->covariance[k*6+l]*work->measurement[j*6+l];
        work->innovation[i*3+j]=value;
    }
    if(tinylinalg_solve(3,work->innovation,residual,solved,work->scratch,72)!=TINYLINALG_OK)return TINYMEKF_ARITHMETIC_ERROR;
    for(size_t i=0;i<3;++i)statistic+=residual[i]*solved[i];
    if(!isfinite(statistic)||statistic<0)return TINYMEKF_ARITHMETIC_ERROR;
    if(gate>0&&statistic>gate)return TINYMEKF_REJECTED;
    for(size_t i=0;i<6;++i){
        double rhs[3]={0};
        for(size_t j=0;j<3;++j)for(size_t k=0;k<6;++k)rhs[j]+=state->covariance[i*6+k]*work->measurement[j*6+k];
        if(tinylinalg_solve(3,work->innovation,rhs,solved,work->scratch,72)!=TINYLINALG_OK)return TINYMEKF_ARITHMETIC_ERROR;
        for(size_t j=0;j<3;++j){work->gain[i*3+j]=solved[j];correction[i]+=solved[j]*residual[j];}
    }
    if(!finite_array(correction,6))return TINYMEKF_ARITHMETIC_ERROR;
    if(tinylinalg_norm3(correction)>.5)return TINYMEKF_REJECTED;
    next=*state;
    for(size_t i=0;i<3;++i){rotation[i]=-correction[i];next.bias[i]+=correction[i+3];if(!isfinite(next.bias[i]))return TINYMEKF_ARITHMETIC_ERROR;}
    if(tinyrotation_from_rotvec(rotation,delta)!=TINYROTATION_OK||tinyrotation_multiply(delta,normalized,next.q)!=TINYROTATION_OK)return TINYMEKF_ARITHMETIC_ERROR;
    identity(work->a);
    for(size_t i=0;i<6;++i)for(size_t j=0;j<6;++j)for(size_t k=0;k<3;++k)work->a[i*6+j]-=work->gain[i*3+k]*work->measurement[k*6+j];
    multiply(work->a,state->covariance,work->b,0);multiply(work->b,work->a,work->c,1);
    for(size_t i=0;i<6;++i)for(size_t j=0;j<6;++j)for(size_t k=0;k<3;++k)work->c[i*6+j]+=variance*work->gain[i*3+k]*work->gain[j*3+k];
    /* Reset the local attitude error after injecting it into the reference. */
    skew(correction,cross);identity(work->a);
    for(size_t i=0;i<3;++i)for(size_t j=0;j<3;++j)work->a[i*6+j]-=.5*cross[i*3+j];
    multiply(work->a,work->c,work->b,0);multiply(work->b,work->a,work->c,1);
    if(!symmetrize(work->c)||!covariance_valid(work->c,work))return TINYMEKF_ARITHMETIC_ERROR;
    memcpy(next.covariance,work->c,sizeof next.covariance);*state=next;if(nis)*nis=statistic;return TINYMEKF_OK;
}
