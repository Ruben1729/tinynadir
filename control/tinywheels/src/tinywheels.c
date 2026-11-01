/* SPDX-License-Identifier: MIT */
#include "tinywheels/tinywheels.h"
#include "tinylinalg/tinylinalg.h"
#include <math.h>
#include <string.h>
tinywheels_status tinywheels_allocate(size_t n,const double *axes,const double requested[3],
    const double *speed,const double *tlimit,const double *slimit,const unsigned char *enabled,
    double *wheel_torque,double residual[3],tinywheels_workspace *work){
    double multiplier[3]={0},out[TINYWHEELS_MAX_WHEELS],error[3];size_t active=0;
    if(n<1||n>TINYWHEELS_MAX_WHEELS||!axes||!requested||!speed||!tlimit||!slimit||!wheel_torque||!residual||!work)return TINYWHEELS_INVALID_ARGUMENT;
    for(size_t i=0;i<3;++i)if(!isfinite(requested[i]))return TINYWHEELS_INVALID_ARGUMENT;
    for(size_t i=0;i<9;++i)work->gram[i]=0;
    for(size_t w=0;w<n;++w){
        if(!isfinite(speed[w])||!isfinite(tlimit[w])||tlimit[w]<0||!isfinite(slimit[w])||slimit[w]<0)return TINYWHEELS_INVALID_ARGUMENT;
        double norm=tinylinalg_norm3(axes+3*w);if(!isfinite(norm)||fabs(norm-1)>1e-10)return TINYWHEELS_INVALID_ARGUMENT;
        if((!enabled||enabled[w])&&tlimit[w]>0&&slimit[w]>0){++active;for(size_t i=0;i<3;++i)for(size_t j=0;j<3;++j)work->gram[i*3+j]+=axes[w*3+i]*axes[w*3+j];}
    }
    if(!active)return TINYWHEELS_NO_AUTHORITY;
    if(tinylinalg_eigen_symmetric(3,work->gram,work->values,work->vectors,work->scratch,18)!=TINYLINALG_OK)return TINYWHEELS_ARITHMETIC_ERROR;
    for(size_t k=0;k<3;++k)if(work->values[k]>work->values[2]*1e-12){
        double projection=0;for(size_t i=0;i<3;++i)projection+=work->vectors[i*3+k]*requested[i];
        projection/=work->values[k];if(!isfinite(projection))return TINYWHEELS_ARITHMETIC_ERROR;
        for(size_t i=0;i<3;++i){multiplier[i]+=work->vectors[i*3+k]*projection;if(!isfinite(multiplier[i]))return TINYWHEELS_ARITHMETIC_ERROR;}
    }
    for(size_t i=0;i<3;++i)error[i]=requested[i];
    for(size_t w=0;w<n;++w){
        out[w]=0;
        if((!enabled||enabled[w])&&tlimit[w]>0&&slimit[w]>0){
            for(size_t i=0;i<3;++i)out[w]-=axes[w*3+i]*multiplier[i];
            if(!isfinite(out[w]))return TINYWHEELS_ARITHMETIC_ERROR;
            out[w]=fmax(-tlimit[w],fmin(tlimit[w],out[w]));
            if((speed[w]>=slimit[w]&&out[w]>0)||(speed[w]<=-slimit[w]&&out[w]<0))out[w]=0;
        }
        for(size_t i=0;i<3;++i)error[i]+=axes[w*3+i]*out[w];
    }
    for(size_t i=0;i<3;++i)if(!isfinite(error[i]))return TINYWHEELS_ARITHMETIC_ERROR;
    memcpy(wheel_torque,out,n*sizeof *wheel_torque);memcpy(residual,error,sizeof error);return TINYWHEELS_OK;
}
