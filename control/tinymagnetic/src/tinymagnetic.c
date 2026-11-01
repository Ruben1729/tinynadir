/* SPDX-License-Identifier: MIT */
#include "tinymagnetic/tinymagnetic.h"
#include "tinylinalg/tinylinalg.h"
#include <math.h>
#include <stddef.h>
#include <string.h>
tinymagnetic_status tinymagnetic_allocate(const double field[3],const double requested[3],
    const double limit[3],double dipole[3],double residual[3]){
    double unit[3],norm,m[3],achieved[3],error[3];
    if(!field||!requested||!limit||!dipole||!residual)return TINYMAGNETIC_INVALID_ARGUMENT;
    for(size_t i=0;i<3;++i)if(!isfinite(field[i])||!isfinite(requested[i])||!isfinite(limit[i])||limit[i]<0)return TINYMAGNETIC_INVALID_ARGUMENT;
    if(tinylinalg_unit3(field,unit)!=TINYLINALG_OK)return TINYMAGNETIC_DEGENERATE;
    norm=tinylinalg_norm3(field);if(!isfinite(norm))return TINYMAGNETIC_ARITHMETIC_ERROR;
    if(tinylinalg_cross3(unit,requested,m)!=TINYLINALG_OK)return TINYMAGNETIC_ARITHMETIC_ERROR;
    for(size_t i=0;i<3;++i){m[i]/=norm;if(!isfinite(m[i]))return TINYMAGNETIC_ARITHMETIC_ERROR;m[i]=fmax(-limit[i],fmin(limit[i],m[i]));}
    if(tinylinalg_cross3(m,field,achieved)!=TINYLINALG_OK)return TINYMAGNETIC_ARITHMETIC_ERROR;
    for(size_t i=0;i<3;++i){error[i]=requested[i]-achieved[i];if(!isfinite(error[i]))return TINYMAGNETIC_ARITHMETIC_ERROR;}
    memcpy(dipole,m,sizeof m);memcpy(residual,error,sizeof error);return TINYMAGNETIC_OK;
}
