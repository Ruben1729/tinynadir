#ifndef TINYRK4_H
#define TINYRK4_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
#define TINYRK4_MAX_STATE 32u
typedef enum { TINYRK4_OK=0,TINYRK4_INVALID,TINYRK4_DERIVATIVE_ERROR,TINYRK4_ARITHMETIC } tinyrk4_status;
/* Callback writes n finite derivatives and returns zero on success. Caller
   controls callback state; library cannot roll back callback side effects. */
typedef int (*tinyrk4_rhs)(double t,const double *state,size_t n,void *context,double *derivative);
/* One fixed RK4 step; state/out may alias. workspace must be disjoint and
   contain >=6*n doubles. Failed calls preserve out; workspace unspecified. */
tinyrk4_status tinyrk4_step(tinyrk4_rhs rhs,void *context,size_t n,double t,double step,const double *state,double *out,double *workspace,size_t count);
#ifdef __cplusplus
}
#endif
#endif
