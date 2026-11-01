/* SPDX-License-Identifier: MIT */
#include "tinyquest/tinyquest.h"
#include "tinylinalg/tinylinalg.h"
#include "tinyrotation/tinyrotation.h"
#include <float.h>
#include <math.h>
#include <string.h>

static tinyquest_status profile(size_t count,const double *reference,const double *body,
    const double *weights,double k[16]){
    double b[9]={0},max_weight=0,total=0,trace;
    if(count<2||count>TINYQUEST_MAX_OBSERVATIONS||!reference||!body)return TINYQUEST_INVALID_ARGUMENT;
    for(size_t i=0;i<count;++i){
        double w=weights?weights[i]:1;
        if(!isfinite(w)||w<0)return TINYQUEST_INVALID_ARGUMENT;
        max_weight=fmax(max_weight,w);
        for(size_t j=0;j<3;++j)if(!isfinite(reference[i*3+j])||!isfinite(body[i*3+j]))return TINYQUEST_INVALID_ARGUMENT;
    }
    if(max_weight==0)return TINYQUEST_DEGENERATE;
    for(size_t i=0;i<count;++i)total+=(weights?weights[i]:1)/max_weight;
    for(size_t i=0;i<count;++i){
        double w=(weights?weights[i]:1)/max_weight/total,r[3],o[3];
        if(w==0)continue;
        if(tinylinalg_unit3(reference+3*i,r)!=TINYLINALG_OK||tinylinalg_unit3(body+3*i,o)!=TINYLINALG_OK)return TINYQUEST_DEGENERATE;
        for(size_t j=0;j<3;++j)for(size_t l=0;l<3;++l)b[j*3+l]+=w*o[j]*r[l];
    }
    trace=b[0]+b[4]+b[8];k[0]=trace;
    k[1]=k[4]=b[7]-b[5];k[2]=k[8]=b[2]-b[6];k[3]=k[12]=b[3]-b[1];
    for(size_t i=0;i<3;++i)for(size_t j=0;j<3;++j)k[(i+1)*4+j+1]=b[i*3+j]+b[j*3+i]-(i==j?trace:0);
    return TINYQUEST_OK;
}
static double minor3(const double m[16],size_t row,size_t col){
    double a[9];size_t n=0;
    for(size_t i=0;i<4;++i)if(i!=row)for(size_t j=0;j<4;++j)if(j!=col)a[n++]=m[i*4+j];
    return a[0]*(a[4]*a[8]-a[5]*a[7])-a[1]*(a[3]*a[8]-a[5]*a[6])+a[2]*(a[3]*a[7]-a[4]*a[6]);
}
static double determinant(const double m[16]){
    double d=0;for(size_t j=0;j<4;++j)d+=(j%2?-1:1)*m[j]*minor3(m,0,j);return d;
}
static int recover(tinyquest_workspace *work,double lambda,double q[4]){
    double best=-1;
    for(size_t i=0;i<4;++i)for(size_t j=0;j<4;++j)work->shifted[i*4+j]=(i==j?lambda:0)-work->profile[i*4+j];
    for(size_t i=0;i<4;++i)for(size_t j=0;j<4;++j)work->adjugate[i*4+j]=((i+j)%2?-1:1)*minor3(work->shifted,j,i);
    for(size_t j=0;j<4;++j){
        double norm=0;for(size_t i=0;i<4;++i)norm+=work->adjugate[i*4+j]*work->adjugate[i*4+j];
        if(norm>best){best=norm;for(size_t i=0;i<4;++i)q[i]=work->adjugate[i*4+j];}
    }
    if(best<1e-24)return 0;
    return tinyrotation_normalize(q,q)==TINYROTATION_OK;
}
tinyquest_status tinyquest_solve(size_t count,const double *reference,
    const double *body,const double *weights,double q[4],double *loss,
    tinyquest_workspace *work){
    double p2=0,p3=0,c2,c3,c4,lambda=1,answer[4],cost,derivative=0;
    int converged=0;tinyquest_status status;
    if(!q||!work)return TINYQUEST_INVALID_ARGUMENT;
    status=profile(count,reference,body,weights,work->profile);if(status!=TINYQUEST_OK)return status;
    /* trace(K)=0. Newton identities give det(lambda*I-K)=
       lambda^4-trace(K^2)/2*lambda^2-trace(K^3)/3*lambda+det(K).
       The largest eigenvalue is <= sum(normalized weights)=1. */
    for(size_t i=0;i<4;++i)for(size_t j=0;j<4;++j){
        p2+=work->profile[i*4+j]*work->profile[j*4+i];
        for(size_t k=0;k<4;++k)p3+=work->profile[i*4+j]*work->profile[j*4+k]*work->profile[k*4+i];
    }
    c2=p2/2;c3=p3/3;c4=determinant(work->profile);
    for(size_t iteration=0;iteration<64;++iteration){
        double value=((lambda*lambda-c2)*lambda-c3)*lambda+c4;
        double step,next;
        derivative=4*lambda*lambda*lambda-2*c2*lambda-c3;
        if(derivative<=1e-12)return TINYQUEST_DEGENERATE;
        step=value/derivative;next=lambda-step;
        if(!isfinite(next)||next<-.001||next>1+1e-12)return TINYQUEST_NO_CONVERGENCE;
        lambda=next;
        if(fabs(step)<=16*DBL_EPSILON){converged=1;break;}
    }
    if(!converged)return TINYQUEST_NO_CONVERGENCE;
    /* Full adjugate columns avoid the Rodrigues/Cayley singularity at 180 deg.
       Rayleigh refinement reduces polynomial cancellation near the optimum. */
    for(size_t refine=0;refine<3;++refine){
        double rayleigh=0;
        if(!recover(work,lambda,answer))return TINYQUEST_DEGENERATE;
        for(size_t i=0;i<4;++i)for(size_t j=0;j<4;++j)rayleigh+=answer[i]*work->profile[i*4+j]*answer[j];
        lambda=rayleigh;
    }
    for(size_t i=0;i<4;++i){
        double residual=-lambda*answer[i];for(size_t j=0;j<4;++j)residual+=work->profile[i*4+j]*answer[j];
        if(fabs(residual)>1e-11)return TINYQUEST_NO_CONVERGENCE;
    }
    if(answer[0]<0)for(size_t i=0;i<4;++i)answer[i]=-answer[i];
    cost=fmax(0,1-lambda);memcpy(q,answer,sizeof answer);if(loss)*loss=cost;
    return TINYQUEST_OK;
}
