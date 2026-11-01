#ifndef TINYGUIDANCE_H
#define TINYGUIDANCE_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYGUIDANCE_OK=0,TINYGUIDANCE_INVALID_ARGUMENT,TINYGUIDANCE_DEGENERATE,TINYGUIDANCE_ARITHMETIC_ERROR} tinyguidance_status;
/* Construct reference->desired-body orientation aligning chosen body pointing
 * axes with two supplied reference directions. Primary has priority; secondary
 * fixes roll. Can describe Sun, inertial, or target pointing. Directions must
 * share an epoch/reference frame and be nonparallel. Failure preserves output. */
tinyguidance_status tinyguidance_point(const double reference_primary[3],const double reference_secondary[3],
    const double body_primary_axis[3],const double body_secondary_axis[3],double q[4]);
/* Inertial r [m], v [m/s]. Desired body z points nadir (-r), x along projected
 * velocity, y completes the right-handed frame (-orbit normal). Returned
 * physical desired body rate is [0,-|r cross v|/|r|²,0]. This instantaneous
 * planar-orbit rate omits plane-precession/perturbation derivatives. */
tinyguidance_status tinyguidance_nadir(const double position[3],const double velocity[3],double q[4],double rate[3]);
/* Target line of sight is target_position-spacecraft_position in one frame;
 * auxiliary reference vector fixes roll. Output is an attitude, not a slew. */
tinyguidance_status tinyguidance_target(const double spacecraft_position[3],const double target_position[3],const double auxiliary_reference[3],
    const double body_primary_axis[3],const double body_secondary_axis[3],double q[4]);
#ifdef __cplusplus
}
#endif
#endif
