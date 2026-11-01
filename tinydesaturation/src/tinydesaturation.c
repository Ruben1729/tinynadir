/* SPDX-License-Identifier: MIT */
#include "tinydesaturation/tinydesaturation.h"
#include "tinymagnetic/tinymagnetic.h"
#include <math.h>
#include <stddef.h>
#include <string.h>
tinydesaturation_status tinydesaturation_command(const double momentum[3],const double target[3],double gain,
    const double field[3],const double limit[3],double dipole[3],double torque[3],double residual[3]){
    double request[3],m[3],r[3],achieved[3];tinymagnetic_status status;
    if(!momentum||!dipole||!torque||!isfinite(gain)||gain<0)return TINYDESATURATION_INVALID_ARGUMENT;
    for(size_t i=0;i<3;++i){if(!isfinite(momentum[i])||(target&&!isfinite(target[i])))return TINYDESATURATION_INVALID_ARGUMENT;request[i]=-gain*(momentum[i]-(target?target[i]:0));if(!isfinite(request[i]))return TINYDESATURATION_ARITHMETIC_ERROR;}
    status=tinymagnetic_allocate(field,request,limit,m,r);
    if(status==TINYMAGNETIC_DEGENERATE)return TINYDESATURATION_DEGENERATE;
    if(status==TINYMAGNETIC_ARITHMETIC_ERROR)return TINYDESATURATION_ARITHMETIC_ERROR;
    if(status!=TINYMAGNETIC_OK)return TINYDESATURATION_INVALID_ARGUMENT;
    for(size_t i=0;i<3;++i){achieved[i]=request[i]-r[i];if(!isfinite(achieved[i]))return TINYDESATURATION_ARITHMETIC_ERROR;}
    memcpy(dipole,m,sizeof m);memcpy(torque,achieved,sizeof achieved);if(residual)memcpy(residual,r,sizeof r);return TINYDESATURATION_OK;
}
