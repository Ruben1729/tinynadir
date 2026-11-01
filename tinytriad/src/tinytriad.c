/* SPDX-License-Identifier: MIT */
#include "tinytriad/tinytriad.h"
#include "tinylinalg/tinylinalg.h"
#include "tinyrotation/tinyrotation.h"
#include <math.h>
#include <string.h>

static tinytriad_status basis(const double first[3],const double second[3],double frame[9]){
    double a[3],b[3],c[3],d[3];tinylinalg_status status;
    status=tinylinalg_unit3(first,a);
    if(status==TINYLINALG_INVALID_ARGUMENT)return TINYTRIAD_INVALID_ARGUMENT;
    if(status!=TINYLINALG_OK)return TINYTRIAD_DEGENERATE;
    status=tinylinalg_unit3(second,b);
    if(status==TINYLINALG_INVALID_ARGUMENT)return TINYTRIAD_INVALID_ARGUMENT;
    if(status!=TINYLINALG_OK)return TINYTRIAD_DEGENERATE;
    if(tinylinalg_cross3(a,b,c)!=TINYLINALG_OK)return TINYTRIAD_ARITHMETIC_ERROR;
    if(tinylinalg_norm3(c)<=1e-10)return TINYTRIAD_DEGENERATE;
    if(tinylinalg_unit3(c,c)!=TINYLINALG_OK)return TINYTRIAD_DEGENERATE;
    if(tinylinalg_cross3(a,c,d)!=TINYLINALG_OK)return TINYTRIAD_ARITHMETIC_ERROR;
    for(size_t i=0;i<3;++i){frame[i*3]=a[i];frame[i*3+1]=c[i];frame[i*3+2]=d[i];}
    return TINYTRIAD_OK;
}
tinytriad_status tinytriad_dcm(const double reference1[3],const double reference2[3],
    const double body1[3],const double body2[3],double dcm[9]){
    double reference[9],body[9],out[9];tinytriad_status status;
    if(!dcm)return TINYTRIAD_INVALID_ARGUMENT;
    status=basis(reference1,reference2,reference);if(status!=TINYTRIAD_OK)return status;
    status=basis(body1,body2,body);if(status!=TINYTRIAD_OK)return status;
    for(size_t i=0;i<3;++i)for(size_t j=0;j<3;++j){out[i*3+j]=0;for(size_t k=0;k<3;++k)out[i*3+j]+=body[i*3+k]*reference[j*3+k];}
    memcpy(dcm,out,sizeof out);return TINYTRIAD_OK;
}
tinytriad_status tinytriad_solve(const double reference1[3],const double reference2[3],
    const double body1[3],const double body2[3],double quaternion[4]){
    double matrix[9];tinytriad_status status;
    if(!quaternion)return TINYTRIAD_INVALID_ARGUMENT;
    status=tinytriad_dcm(reference1,reference2,body1,body2,matrix);if(status!=TINYTRIAD_OK)return status;
    if(tinyrotation_from_dcm(matrix,quaternion)!=TINYROTATION_OK)return TINYTRIAD_ARITHMETIC_ERROR;
    return TINYTRIAD_OK;
}
