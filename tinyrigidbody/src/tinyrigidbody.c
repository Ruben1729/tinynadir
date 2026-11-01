/* SPDX-License-Identifier: MIT */
#include "tinyrigidbody/tinyrigidbody.h"
#include "tinylinalg/tinylinalg.h"
#include "tinygyro/tinygyro.h"
#include <math.h>
#include <stddef.h>
#include <string.h>
tinyrigidbody_status tinyrigidbody_derivative(const double q[4],const double omega[3],const double inertia[9],
    const double torque[3],const double wheels[3],double qdot[4],double omega_dot[3],tinyrigidbody_workspace *work){
    double momentum[3],cross[3],rhs[3],dq[4],dw[3];tinylinalg_status status;
    if(!q||!omega||!inertia||!torque||!qdot||!omega_dot||!work)return TINYRIGIDBODY_INVALID_ARGUMENT;
    for(size_t i=0;i<3;++i)if(!isfinite(omega[i])||!isfinite(torque[i])||(wheels&&!isfinite(wheels[i])))return TINYRIGIDBODY_INVALID_ARGUMENT;
    status=tinylinalg_cholesky(3,inertia,work->factor,work->scratch,12);
    if(status!=TINYLINALG_OK)return TINYRIGIDBODY_INVALID_INERTIA;
    for(size_t i=0;i<3;++i){momentum[i]=wheels?wheels[i]:0;for(size_t j=0;j<3;++j)momentum[i]+=inertia[i*3+j]*omega[j];if(!isfinite(momentum[i]))return TINYRIGIDBODY_ARITHMETIC_ERROR;}
    if(tinylinalg_cross3(omega,momentum,cross)!=TINYLINALG_OK)return TINYRIGIDBODY_ARITHMETIC_ERROR;
    for(size_t i=0;i<3;++i){rhs[i]=torque[i]-cross[i];if(!isfinite(rhs[i]))return TINYRIGIDBODY_ARITHMETIC_ERROR;}
    if(tinylinalg_solve(3,inertia,rhs,dw,work->scratch,12)!=TINYLINALG_OK)return TINYRIGIDBODY_ARITHMETIC_ERROR;
    if(tinygyro_derivative(q,omega,NULL,dq)!=TINYGYRO_OK)return TINYRIGIDBODY_INVALID_ARGUMENT;
    memcpy(qdot,dq,sizeof dq);memcpy(omega_dot,dw,sizeof dw);return TINYRIGIDBODY_OK;
}
