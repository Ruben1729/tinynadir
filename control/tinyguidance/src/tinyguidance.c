/* SPDX-License-Identifier: MIT */
#include "tinyguidance/tinyguidance.h"
#include "tinytriad/tinytriad.h"
#include "tinylinalg/tinylinalg.h"
#include <math.h>
#include <stddef.h>
#include <string.h>
tinyguidance_status tinyguidance_point(const double reference_primary[3],const double reference_secondary[3],
    const double body_primary[3],const double body_secondary[3],double q[4]){
    tinytriad_status status=tinytriad_solve(reference_primary,reference_secondary,body_primary,body_secondary,q);
    if(status==TINYTRIAD_INVALID_ARGUMENT)return TINYGUIDANCE_INVALID_ARGUMENT;
    if(status==TINYTRIAD_DEGENERATE)return TINYGUIDANCE_DEGENERATE;
    if(status!=TINYTRIAD_OK)return TINYGUIDANCE_ARITHMETIC_ERROR;
    return TINYGUIDANCE_OK;
}
tinyguidance_status tinyguidance_nadir(const double position[3],const double velocity[3],double q[4],double rate[3]){
    const double body_primary[3]={0,0,1},body_secondary[3]={1,0,0};
    double nadir[3],h[3],out[4],w[3]={0},radius;
    if(!q||!rate||!position||!velocity)return TINYGUIDANCE_INVALID_ARGUMENT;
    for(size_t i=0;i<3;++i)if(!isfinite(position[i])||!isfinite(velocity[i]))return TINYGUIDANCE_INVALID_ARGUMENT;
    if(tinylinalg_unit3(position,nadir)!=TINYLINALG_OK)return TINYGUIDANCE_DEGENERATE;
    if(tinylinalg_cross3(nadir,velocity,h)!=TINYLINALG_OK)return TINYGUIDANCE_ARITHMETIC_ERROR;
    radius=tinylinalg_norm3(position);if(!isfinite(radius))return TINYGUIDANCE_ARITHMETIC_ERROR;
    w[1]=-tinylinalg_norm3(h)/radius;if(!isfinite(w[1]))return TINYGUIDANCE_ARITHMETIC_ERROR;
    for(size_t i=0;i<3;++i)nadir[i]=-nadir[i];
    tinyguidance_status status=tinyguidance_point(nadir,velocity,body_primary,body_secondary,out);
    if(status!=TINYGUIDANCE_OK)return status;
    memcpy(q,out,sizeof out);memcpy(rate,w,sizeof w);return TINYGUIDANCE_OK;
}
tinyguidance_status tinyguidance_target(const double spacecraft[3],const double target[3],const double auxiliary[3],
    const double body_primary[3],const double body_secondary[3],double q[4]){
    double line[3];
    if(!spacecraft||!target)return TINYGUIDANCE_INVALID_ARGUMENT;
    for(size_t i=0;i<3;++i){if(!isfinite(spacecraft[i])||!isfinite(target[i]))return TINYGUIDANCE_INVALID_ARGUMENT;line[i]=target[i]-spacecraft[i];if(!isfinite(line[i]))return TINYGUIDANCE_ARITHMETIC_ERROR;}
    return tinyguidance_point(line,auxiliary,body_primary,body_secondary,q);
}
