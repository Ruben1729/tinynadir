/* SPDX-License-Identifier: MIT */
#ifndef TINYCOVARIANCE_H
#define TINYCOVARIANCE_H
#include <stddef.h>
#include "tinyelements/tinyelements.h"
#ifdef __cplusplus
extern "C" {
#endif
#define TINYCOVARIANCE_MAX_DIM 6
typedef enum { TINYCOVARIANCE_OK,TINYCOVARIANCE_INVALID,TINYCOVARIANCE_NOT_PSD,TINYCOVARIANCE_MODEL_ERROR,TINYCOVARIANCE_SINGULAR } tinycovariance_status;
typedef int (*tinycovariance_model)(const double *input,size_t input_dimension,double *output,size_t output_dimension,void *context);
/* Matrices row-major, dimensions 1..6. P must be finite symmetric PSD.
   Computes Pout=J*P*J^T. All failures preserve output. */
tinycovariance_status tinycovariance_transform(size_t input_dimension,size_t output_dimension,const double *jacobian,const double *covariance,double *output);
/* Central finite differences. Each positive step is in the input's physical
   unit. output_period may be NULL or contain zero/positive periods (e.g. 2*pi)
   to unwrap angular differences. Callback writes all outputs, returns nonzero.
   At most 2*input_dimension calls; no allocation. */
tinycovariance_status tinycovariance_jacobian(size_t input_dimension,size_t output_dimension,
 const double *input,const double *step,const double *output_period,tinycovariance_model model,void *context,double *jacobian);
/* A fixed known rotation applied to both position and velocity coordinates.
   Does not differentiate an orbit-dependent frame or include its angular rate. */
tinycovariance_status tinycovariance_rotate_state(const double rotation[9],const double covariance[36],double output[36]);
/* Convenience orbital maps use central finite-difference Jacobians.
   Classical vector ordered [p_m,e,i,Omega,omega,nu], not semimajor axis/mean anomaly.
   Classical maps reject e<1e-8 or |sin(i)|<1e-8. Equinoctial convention is
   tinyelements' prograde modified [p_m,f,g,h,k,L] convention. SI/radians. */
tinycovariance_status tinycovariance_cartesian_to_classical(double mu,const tinyelements_state *state,const double cartesian_covariance[36],tinyelements_classical *elements,double classical_covariance[36]);
tinycovariance_status tinycovariance_classical_to_cartesian(double mu,const tinyelements_classical *elements,const double classical_covariance[36],tinyelements_state *state,double cartesian_covariance[36]);
tinycovariance_status tinycovariance_cartesian_to_equinoctial(double mu,const tinyelements_state *state,const double cartesian_covariance[36],tinyelements_equinoctial *elements,double equinoctial_covariance[36]);
tinycovariance_status tinycovariance_equinoctial_to_cartesian(double mu,const tinyelements_equinoctial *elements,const double equinoctial_covariance[36],tinyelements_state *state,double cartesian_covariance[36]);
#ifdef __cplusplus
}
#endif
#endif
