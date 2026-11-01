#ifndef TINYMAGNETIC_H
#define TINYMAGNETIC_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYMAGNETIC_OK=0,TINYMAGNETIC_INVALID_ARGUMENT,TINYMAGNETIC_DEGENERATE,TINYMAGNETIC_ARITHMETIC_ERROR} tinymagnetic_status;
/* Minimum-norm dipole m=(B cross requested_torque)/|B|^2 followed by per-axis
 * clipping. B is body magnetic field [T], torque [N m], limits [A m^2].
 * Returns dipole and residual=requested-cross(dipole,B). Torque parallel to B
 * is unachievable. Outputs must be distinct and are preserved on failure. */
tinymagnetic_status tinymagnetic_allocate(const double field[3],const double requested_torque[3],
    const double dipole_limit[3],double dipole[3],double residual[3]);
#ifdef __cplusplus
}
#endif
#endif
