/* SPDX-License-Identifier: MIT */
#ifndef TINYBATCHLS_H
#define TINYBATCHLS_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
#define TINYBATCHLS_MAX_PARAMETERS 8
#define TINYBATCHLS_MAX_OBSERVATIONS 64
typedef enum {TINYBATCHLS_OK,TINYBATCHLS_INVALID,TINYBATCHLS_MODEL_ERROR,TINYBATCHLS_RANK_DEFICIENT,TINYBATCHLS_NO_CONVERGENCE} tinybatchls_status;
typedef int (*tinybatchls_model)(const double *parameters,size_t parameter_count,double *predictions,size_t observation_count,void *context);
typedef struct {
    size_t parameter_count,observation_count,max_iterations;
    double parameter_scale[8],relative_difference_step,step_tolerance,gradient_tolerance,cost_tolerance,initial_damping;
} tinybatchls_options;
typedef struct {
    double parameters[8],covariance[64],weighted_sum_squares;
    size_t iterations,model_evaluations;
} tinybatchls_result;
/* Caller-owned scratch; contents unspecified after each call, including failure. */
typedef struct {
    double jacobian[64*8],prediction[64],trial_prediction[64],plus[64],minus[64];
    double normal[64],factor[64],gradient[8],delta[8],trial[8];
} tinybatchls_workspace;
tinybatchls_status tinybatchls_default_options(size_t parameter_count,size_t observation_count,tinybatchls_options *options);
/* Weighted Levenberg-Marquardt with central numerical Jacobians.
   sigma gives strictly positive known measurement standard deviations;
   covariance=(J^T*W*J)^-1 in physical parameter units, not residual-scaled.
   parameter_scale defines positive physical scaling. Callback returns nonzero
   and writes every prediction. At most 100 outer iterations, 8 parameters,
   64 scalar observations. No parameter bounds or robust loss is implicit.
   Scratch must not overlap inputs/output. Failure preserves result. */
tinybatchls_status tinybatchls_solve(const tinybatchls_options *options,
 tinybatchls_model model,void *context,const double *observations,const double *sigma,
 const double *initial_parameters,tinybatchls_workspace *workspace,tinybatchls_result *result);
#ifdef __cplusplus
}
#endif
#endif
