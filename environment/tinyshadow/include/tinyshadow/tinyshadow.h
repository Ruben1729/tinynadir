#ifndef TINYSHADOW_H
#define TINYSHADOW_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYSHADOW_OK=0,TINYSHADOW_INVALID,TINYSHADOW_GEOMETRY,TINYSHADOW_ARITHMETIC,TINYSHADOW_UNBRACKETED,TINYSHADOW_NO_CONVERGENCE,TINYSHADOW_CALLBACK_ERROR} tinyshadow_status;
/* Spherical bodies, apparent angular-circle overlap profile; SI positions
   relative to occulting body's center. No refraction or atmospheric shadow.
   Output fraction0..1, failed calls preserve outputs. */
tinyshadow_status tinyshadow_sun_fraction(const double spacecraft_m[3],const double sun_m[3],double body_radius_m,double sun_radius_m,double *fraction);
/* Parallel sunlight/cylindrical shadow approximation. */
tinyshadow_status tinyshadow_cylindrical(const double spacecraft_m[3],const double sun_direction[3],double body_radius_m,int *eclipsed);
/* Axis-aligned ellipsoid at origin. Tangent or inside contacts block LOS. */
tinyshadow_status tinyshadow_segment_ellipsoid(const double from_m[3],const double to_m[3],const double semiaxes_m[3],int *blocked);
typedef int (*tinyshadow_fraction_callback)(double time_s,void *context,double *fraction);
/* Find ONE bracketed threshold crossing. Does not search for or assume that
   all crossings were found. At most max_iterations+2 callback evaluations. */
tinyshadow_status tinyshadow_event(tinyshadow_fraction_callback callback,void *context,double start_s,double end_s,double threshold,double tolerance_s,unsigned max_iterations,double *event_s);
#ifdef __cplusplus
}
#endif
#endif
