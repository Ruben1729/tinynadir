#ifndef TINYINTERPOLATION_H
#define TINYINTERPOLATION_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef enum { TINYINTERPOLATION_OK=0, TINYINTERPOLATION_INVALID, TINYINTERPOLATION_RANGE, TINYINTERPOLATION_ARITHMETIC } tinyinterpolation_status;
/* Strictly increasing finite abscissae; no extrapolation. All failed calls
   preserve outputs. Linear tables: 2..4096 points; Lagrange: 2..8 points. */
tinyinterpolation_status tinyinterpolation_linear(const double *x,const double *y,size_t n,double query,double *value);
tinyinterpolation_status tinyinterpolation_lagrange(const double *x,const double *y,size_t n,double query,double *value);
/* Endpoint slopes use dy/dx. derivative may be NULL; other outputs must be
   distinct. Evaluation only within [x0,x1]. */
tinyinterpolation_status tinyinterpolation_hermite(double x0,double x1,double y0,double y1,double slope0,double slope1,double query,double *value,double *derivative);
#ifdef __cplusplus
}
#endif
#endif
