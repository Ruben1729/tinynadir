#ifndef TINYELLIPTIC_H
#define TINYELLIPTIC_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYELLIPTIC_OK=0,TINYELLIPTIC_INVALID_ARGUMENT,TINYELLIPTIC_RANGE,TINYELLIPTIC_SINGULAR,TINYELLIPTIC_NO_CONVERGENCE,TINYELLIPTIC_ARITHMETIC_ERROR} tinyelliptic_status;
/* Real Legendre integrals with parameter m=k² in [0,1], amplitudes in radians,
 * |amplitude| <=1e6. These use parameter m, not elliptic modulus k.
 * First integral is singular at m=1 across |amplitude|>=pi/2. Outputs atomic. */
tinyelliptic_status tinyelliptic_first(double amplitude,double m,double *value);
tinyelliptic_status tinyelliptic_second(double amplitude,double m,double *value);
tinyelliptic_status tinyelliptic_complete_first(double m,double *value);
tinyelliptic_status tinyelliptic_complete_second(double m,double *value);
/* Invert the monotonic real second integral. Result amplitude must be <=1e6
 * in magnitude; ill conditioning near m=1, amplitude=pi/2 is intrinsic. */
tinyelliptic_status tinyelliptic_inverse_second(double value,double m,double *amplitude);
/* Signed ellipse arc for x=a*cos(theta),y=b*sin(theta), a,b>0 in any common
 * length unit. Endpoint amplitudes bounded by 1e6 radians. */
tinyelliptic_status tinyelliptic_arc(double a,double b,double start,double end,double *length);
#ifdef __cplusplus
}
#endif
#endif
