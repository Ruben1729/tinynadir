#ifndef TINYDESATURATION_H
#define TINYDESATURATION_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYDESATURATION_OK=0,TINYDESATURATION_INVALID_ARGUMENT,TINYDESATURATION_DEGENERATE,TINYDESATURATION_ARITHMETIC_ERROR} tinydesaturation_status;
/* Request external body torque -gain*(momentum-target) [N m], then perform
 * magnetic allocation. momentum/optional target are body-coordinate wheel
 * angular momentum [N m s], gain>=0 [1/s], field [T], limits [A m^2].
 * Returns dipole and actual external unloading torque; optional residual is
 * unachievable requested torque. Output buffers must be distinct. The caller
 * coordinates reaction-wheel counter-torque with the attitude controller. */
tinydesaturation_status tinydesaturation_command(const double momentum[3],const double target[3],double gain,
    const double field[3],const double dipole_limit[3],double dipole[3],double torque[3],double residual[3]);
#ifdef __cplusplus
}
#endif
#endif
