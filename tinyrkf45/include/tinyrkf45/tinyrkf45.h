#ifndef TINYRKF45_H
#define TINYRKF45_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
#define TINYRKF45_MAX_STATE 32u
typedef enum {TINYRKF45_OK=0,TINYRKF45_INVALID,TINYRKF45_DERIVATIVE_ERROR,TINYRKF45_STEP_LIMIT,TINYRKF45_STEP_UNDERFLOW,TINYRKF45_ARITHMETIC} tinyrkf45_status;
typedef int (*tinyrkf45_rhs)(double t,const double *state,size_t n,void *context,double *derivative);
typedef struct {double absolute_tolerance,relative_tolerance,initial_step,min_step,max_step;size_t max_attempts;} tinyrkf45_options;
typedef struct {size_t accepted,rejected,evaluations;double final_time,last_step;} tinyrkf45_stats;
/* Adaptive classical Fehlberg 4(5), accepting the fifth-order state. 1..32
   components; workspace >=9*n doubles, disjoint from state/out/stats. State
   and out may alias. All outputs preserved on failure; scratch and callback
   side effects are not rolled back. max_attempts <=1000000. */
tinyrkf45_status tinyrkf45_integrate(tinyrkf45_rhs rhs,void *context,size_t n,double start,double end,const double *state,const tinyrkf45_options *options,double *out,tinyrkf45_stats *stats,double *workspace,size_t count);
#ifdef __cplusplus
}
#endif
#endif
