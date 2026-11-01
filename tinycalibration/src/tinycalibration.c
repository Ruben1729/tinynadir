/* SPDX-License-Identifier: MIT */
#include "tinycalibration/tinycalibration.h"
#include "tinyrotation/tinyrotation.h"
#include <math.h>
#include <stddef.h>
#include <string.h>
tinycalibration_status tinycalibration_apply(const tinycalibration_parameters *p,const double raw[3],double body[3]){
    double centered[3],sensor[3],out[3],q[4];
    if(!p||!raw||!body)return TINYCALIBRATION_INVALID_ARGUMENT;
    for(size_t i=0;i<3;++i){if(!isfinite(raw[i])||!isfinite(p->bias[i]))return TINYCALIBRATION_INVALID_ARGUMENT;centered[i]=raw[i]-p->bias[i];}
    for(size_t i=0;i<9;++i)if(!isfinite(p->correction[i]))return TINYCALIBRATION_INVALID_ARGUMENT;
    if(tinyrotation_from_dcm(p->sensor_to_body,q)!=TINYROTATION_OK)return TINYCALIBRATION_INVALID_ALIGNMENT;
    for(size_t i=0;i<3;++i){sensor[i]=0;for(size_t j=0;j<3;++j)sensor[i]+=p->correction[3*i+j]*centered[j];}
    for(size_t i=0;i<3;++i){out[i]=0;for(size_t j=0;j<3;++j)out[i]+=p->sensor_to_body[3*i+j]*sensor[j];if(!isfinite(out[i]))return TINYCALIBRATION_ARITHMETIC_ERROR;}
    memcpy(body,out,sizeof out);return TINYCALIBRATION_OK;
}
