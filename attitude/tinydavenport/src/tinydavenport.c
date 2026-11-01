/* SPDX-License-Identifier: MIT */
#include "tinydavenport/tinydavenport.h"
#include "tinylinalg/tinylinalg.h"
#include "tinyrotation/tinyrotation.h"
#include <math.h>
#include <string.h>

static tinydavenport_status profile(size_t count,const double *reference,const double *body,
    const double *weights,double k[16]){
    double b[9]={0},max_weight=0,total=0,trace;
    if(count<2||count>TINYDAVENPORT_MAX_OBSERVATIONS||!reference||!body)return TINYDAVENPORT_INVALID_ARGUMENT;
    for(size_t i=0;i<count;++i){
        double w=weights?weights[i]:1;
        if(!isfinite(w)||w<0)return TINYDAVENPORT_INVALID_ARGUMENT;
        max_weight=fmax(max_weight,w);
        for(size_t j=0;j<3;++j)if(!isfinite(reference[i*3+j])||!isfinite(body[i*3+j]))return TINYDAVENPORT_INVALID_ARGUMENT;
    }
    if(max_weight==0)return TINYDAVENPORT_DEGENERATE;
    for(size_t i=0;i<count;++i)total+=(weights?weights[i]:1)/max_weight;
    for(size_t i=0;i<count;++i){
        double w=(weights?weights[i]:1)/max_weight/total,r[3],o[3];
        if(w==0)continue;
        if(tinylinalg_unit3(reference+3*i,r)!=TINYLINALG_OK||tinylinalg_unit3(body+3*i,o)!=TINYLINALG_OK)return TINYDAVENPORT_DEGENERATE;
        for(size_t j=0;j<3;++j)for(size_t l=0;l<3;++l)b[j*3+l]+=w*o[j]*r[l];
    }
    trace=b[0]+b[4]+b[8];k[0]=trace;
    /* Hamilton/reference-to-body convention: Z is minus the frequently
       published Shuster/passive-convention vector. */
    k[1]=k[4]=b[7]-b[5];k[2]=k[8]=b[2]-b[6];k[3]=k[12]=b[3]-b[1];
    for(size_t i=0;i<3;++i)for(size_t j=0;j<3;++j)k[(i+1)*4+j+1]=b[i*3+j]+b[j*3+i]-(i==j?trace:0);
    return TINYDAVENPORT_OK;
}
tinydavenport_status tinydavenport_solve(size_t count,const double *reference,
    const double *body,const double *weights,double q[4],double *loss,
    tinydavenport_workspace *work){
    double answer[4],cost;tinydavenport_status status;
    if(!q||!work)return TINYDAVENPORT_INVALID_ARGUMENT;
    status=profile(count,reference,body,weights,work->profile);if(status!=TINYDAVENPORT_OK)return status;
    if(tinylinalg_eigen_symmetric(4,work->profile,work->eigenvalues,work->eigenvectors,work->scratch,32)!=TINYLINALG_OK)return TINYDAVENPORT_NO_CONVERGENCE;
    if(work->eigenvalues[3]-work->eigenvalues[2]<=1e-12)return TINYDAVENPORT_DEGENERATE;
    for(size_t i=0;i<4;++i)answer[i]=work->eigenvectors[i*4+3];
    if(answer[0]<0)for(size_t i=0;i<4;++i)answer[i]=-answer[i];
    if(tinyrotation_normalize(answer,answer)!=TINYROTATION_OK)return TINYDAVENPORT_NO_CONVERGENCE;
    cost=fmax(0,1-work->eigenvalues[3]);
    memcpy(q,answer,sizeof answer);if(loss)*loss=cost;return TINYDAVENPORT_OK;
}
