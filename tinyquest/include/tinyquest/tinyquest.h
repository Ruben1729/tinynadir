#ifndef TINYQUEST_H
#define TINYQUEST_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
#define TINYQUEST_MAX_OBSERVATIONS 64u
typedef enum {TINYQUEST_OK=0,TINYQUEST_INVALID_ARGUMENT,TINYQUEST_DEGENERATE,TINYQUEST_NO_CONVERGENCE} tinyquest_status;
typedef struct {double profile[16],shifted[16],adjugate[16];} tinyquest_workspace;
/* QUEST weighted Wahba solution. Arrays contain count consecutive [x,y,z]
 * reference/body directions, count in [2,64]. NULL weights select equal weights;
 * otherwise finite nonnegative weights are rescaled and normalized internally.
 * Positive-weight directions must be nonzero; all components must be finite.
 * q is scalar-first Hamilton, body=C(q)*reference. Optional loss is normalized
 * Wahba loss. Workspace must be separate from all inputs/outputs. Failed calls
 * preserve q/loss. At most 64 quartic Newton steps and three adjugate refinements.
 * Degenerate/poorly separated optimal solutions return a status, not an attitude. */
tinyquest_status tinyquest_solve(size_t count,const double *reference,
    const double *body,const double *weights,double q[4],double *loss,
    tinyquest_workspace *workspace);
#ifdef __cplusplus
}
#endif
#endif
