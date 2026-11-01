#ifndef TINYDAVENPORT_H
#define TINYDAVENPORT_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
#define TINYDAVENPORT_MAX_OBSERVATIONS 64u
typedef enum {TINYDAVENPORT_OK=0,TINYDAVENPORT_INVALID_ARGUMENT,TINYDAVENPORT_DEGENERATE,TINYDAVENPORT_NO_CONVERGENCE} tinydavenport_status;
typedef struct {double profile[16],eigenvalues[4],eigenvectors[16],scratch[32];} tinydavenport_workspace;
/* Flattened arrays contain count consecutive [x,y,z] direction vectors.
 * count is in [2,64]. weights may be NULL for equal weights; otherwise weights
 * are finite, nonnegative and need not sum to one. Zero-weight pairs are ignored.
 * Positive-weight vectors must be nonzero; all components must be finite.
 * q is scalar-first Hamilton, body=C(q)*reference. Optional loss is normalized
 * Wahba loss 0.5*sum(w_i*|b_i-C*r_i|^2)/sum(w_i).
 * Workspace must not overlap inputs/outputs. q/loss are preserved on failure;
 * scratch contents are unspecified. A normalized eigenvalue gap <=1e-12 is
 * rejected as numerically unobservable/nonunique. */
tinydavenport_status tinydavenport_solve(size_t count,const double *reference,
    const double *body,const double *weights,double q[4],double *loss,
    tinydavenport_workspace *workspace);
#ifdef __cplusplus
}
#endif
#endif
