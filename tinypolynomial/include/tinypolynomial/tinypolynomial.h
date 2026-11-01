#ifndef TINYPOLYNOMIAL_H
#define TINYPOLYNOMIAL_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
#define TINYPOLYNOMIAL_MAX_DEGREE 8u
typedef struct {double real,imaginary;} tinypolynomial_complex;
typedef enum {TINYPOLYNOMIAL_OK=0,TINYPOLYNOMIAL_INVALID,TINYPOLYNOMIAL_NO_CONVERGENCE,TINYPOLYNOMIAL_ARITHMETIC} tinypolynomial_status;
/* Real coefficients, descending powers, nonzero leading coefficient.
   Degree1..8. Workspace >=2*degree complex entries, disjoint from all other
   arrays. All roots returned sorted lexicographically; conjugate ordering
   should not be used to identify physical roots. Outputs preserved on failure.
   Repeated/clustered roots are ill-conditioned; small backward residual is
   not a forward-error guarantee. Maximum256 Aberth iterations. */
tinypolynomial_status tinypolynomial_roots(unsigned degree,const double *coefficients,tinypolynomial_complex *roots,tinypolynomial_complex *workspace,size_t workspace_count);
#ifdef __cplusplus
}
#endif
#endif
