#ifndef TINYGRAVITY_H
#define TINYGRAVITY_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
#define TINYGRAVITY_MAX_DEGREE 32u
typedef enum {TINYGRAVITY_OK=0,TINYGRAVITY_INVALID,TINYGRAVITY_SINGULAR,TINYGRAVITY_ARITHMETIC} tinygravity_status;
/* Packed triangular index n*(n+1)/2+m. Coefficients are unnormalized,
   no Condon-Shortley phase; C00=1 includes central gravity. SI body-fixed
   position and acceleration; potential is positive mu/r. No field data is
   implicitly installed. Caller chooses coefficient set, degree and order. */
size_t tinygravity_coefficient_count(unsigned degree);
/* Workspace is 8*coefficient_count doubles, disjoint from inputs/outputs. */
size_t tinygravity_workspace_count(unsigned degree);
tinygravity_status tinygravity_evaluate(double mu_m3_s2,double reference_radius_m,unsigned degree,unsigned order,const double *cosine,const double *sine,size_t coefficient_count,const double position_m[3],double acceleration_m_s2[3],double *potential_m2_s2,double *workspace,size_t workspace_count);
/* Fully normalized geodesy coefficients -> unnormalized coefficients.
   Input/output arrays may alias; outputs preserved on invalid input. */
tinygravity_status tinygravity_unnormalize(unsigned degree,const double *normalized_c,const double *normalized_s,size_t count,double *c,double *s);
#ifdef __cplusplus
}
#endif
#endif
