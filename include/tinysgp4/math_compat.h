#ifndef TGP4_MATH_COMPAT_H
#define TGP4_MATH_COMPAT_H
/* A freestanding firmware build need not provide hosted C library headers.
   These standard C math functions must still be provided by its math library.
   Hosted builds use the platform's normal declarations and isfinite macro. */
#if defined(__STDC_HOSTED__) && !__STDC_HOSTED__
double sin(double);
double cos(double);
double sqrt(double);
double pow(double, double);
double floor(double);
double fabs(double);
double fmod(double, double);
double atan2(double, double);
#define TGP4_ISFINITE(x) __builtin_isfinite(x)
#else
#include <math.h>
#define TGP4_ISFINITE(x) isfinite(x)
#endif
#endif
