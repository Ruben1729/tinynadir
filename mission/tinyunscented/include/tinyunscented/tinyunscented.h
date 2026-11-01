/* SPDX-License-Identifier: MIT */
#ifndef TINYUNSCENTED_H
#define TINYUNSCENTED_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
#define TINYUNSCENTED_MAX_DIM 6
#define TINYUNSCENTED_MAX_POINTS 13
typedef enum { TINYUNSCENTED_OK,TINYUNSCENTED_INVALID,TINYUNSCENTED_NOT_PSD,TINYUNSCENTED_MODEL_ERROR } tinyunscented_status;
typedef struct {
    size_t dimension,count;
    /* Packed point rows: point k starts at points[k*dimension]. */
    double points[TINYUNSCENTED_MAX_POINTS*TINYUNSCENTED_MAX_DIM];
    double mean_weights[TINYUNSCENTED_MAX_POINTS],covariance_weights[TINYUNSCENTED_MAX_POINTS];
} tinyunscented_sigma;
typedef int (*tinyunscented_model)(const double *input,size_t input_dimension,
                                 double *output,size_t output_dimension,void *context);
/* Scaled Merwe sigma points, alpha>0, beta>=0, n+kappa>0. Covariance
   is row-major n*n, symmetric PSD. Point ordering: mean, +columns(L),
   -columns(L), with L*L^T=alpha^2*(n+kappa)*P. All failures preserve outputs. */
tinyunscented_status tinyunscented_generate(size_t n,const double *mean,const double *covariance,
                                          double alpha,double beta,double kappa,tinyunscented_sigma *output);
/* Weighted Euclidean statistics; points packed row-major. Signed weights are
   allowed and mean weights must sum to one. No angle/quaternion/manifold
   averaging or process noise is implicit. Arbitrary signed weights need not
   produce PSD results. */
tinyunscented_status tinyunscented_statistics(size_t count,size_t dimension,const double *points,
 const double *mean_weights,const double *covariance_weights,double *mean,double *covariance);
/* Callback returns nonzero on success and writes every output element.
   At most 13 calls; no hidden allocation. Output dimension 1..6. */
tinyunscented_status tinyunscented_transform(const tinyunscented_sigma *sigma,size_t output_dimension,
 tinyunscented_model model,void *context,double *mean,double *covariance);
#ifdef __cplusplus
}
#endif
#endif
