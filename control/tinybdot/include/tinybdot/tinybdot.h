#ifndef TINYBDOT_H
#define TINYBDOT_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYBDOT_OK=0,TINYBDOT_INVALID_ARGUMENT,TINYBDOT_ARITHMETIC_ERROR} tinybdot_status;
typedef struct {double previous_field[3],filtered_derivative[3];int initialized;} tinybdot_state;
/* Zero-initialize state. field is calibrated body-frame B [T], dt positive [s].
 * Gain >=0 maps T/s to A m^2; time_constant>=0 gives a backward-Euler low pass
 * alpha=dt/(time_constant+dt). Commands are clipped to per-axis dipole limits.
 * First sample initializes state and commands zero. limited_axes may be NULL.
 * Failures preserve state/outputs. Sample with torquers suitably de-energized. */
tinybdot_status tinybdot_step(tinybdot_state *state,const double field[3],double dt,
    double gain,double time_constant,const double dipole_limit[3],double dipole[3],unsigned *limited_axes);
#ifdef __cplusplus
}
#endif
#endif
