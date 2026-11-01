#ifndef TINYDOUBLER_H
#define TINYDOUBLER_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYDOUBLER_OK=0,TINYDOUBLER_INVALID_ARGUMENT,TINYDOUBLER_GEOMETRY,TINYDOUBLER_NO_CONVERGENCE} tinydoubler_status;
typedef struct {double position_m[3],velocity_m_s[3],ranges_m[3],radii_m[2],time_residual_s[2];unsigned iterations;} tinydoubler_result;
/* Distinct double-r method: solve two geocentric radii using Newton updates of
 * conic time-of-flight residuals; ranges1/2 use the positive far intersection
 * of each ray/radius sphere and range3 enforces coplanarity.
 * Elliptic profile e<0.95, nonradial minor arcs, each TOF below half a period.
 * Observations are three consecutive inertial LOS/site xyz triples; SI units.
 * tau1<0<tau3 are relative seconds, magnitudes<=1e6. Positive initial geocentric
 * radii must yield an admissible elliptic geometry. time_tolerance in [1e-10,1]
 * seconds; max_iterations in [1,64]. Up to12 line-search probes/iteration.
 * Multiple solutions depend on seeds; no automatic mission prior is imposed.
 * No heap or persistent scratch; failure preserves result. */
tinydoubler_status tinydoubler_estimate(double mu,double tau1,double tau3,const double line_of_sight[9],
    const double observer_position_m[9],double initial_radius1_m,double initial_radius2_m,
    double time_tolerance_s,unsigned max_iterations,tinydoubler_result *result);
#ifdef __cplusplus
}
#endif
#endif
