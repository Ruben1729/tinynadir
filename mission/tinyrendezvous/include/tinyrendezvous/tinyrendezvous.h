#ifndef TINYRENDEZVOUS_H
#define TINYRENDEZVOUS_H
#include "tinykepler/tinykepler.h"
#include "tinylambert/tinylambert.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYRENDEZVOUS_OK=0,TINYRENDEZVOUS_INVALID,TINYRENDEZVOUS_RANGE,TINYRENDEZVOUS_GEOMETRY,TINYRENDEZVOUS_PROPAGATION,TINYRENDEZVOUS_LAMBERT} tinyrendezvous_status;
typedef struct {double duration_s,phasing_period_s,semimajor_axis_m,other_apsis_radius_m,departure_tangent_delta_v_m_s,arrival_tangent_delta_v_m_s,total_delta_v_m_s;} tinyrendezvous_phasing;
/* Same circular coplanar orbit, target phase=target−chaser in [-pi,+pi].
   Rendezvous at the departure point: target ends at angle2*pi*target_turns,
   chaser completes chaser_turns ellipses. minimum_radius_m is an explicit
   central-body/altitude exclusion radius (0 disables that physical clearance). */
tinyrendezvous_status tinyrendezvous_phase_circular(double mu_m3_s2,double radius_m,double target_phase_rad,unsigned target_turns,unsigned chaser_turns,double minimum_radius_m,tinyrendezvous_phasing *out);
typedef struct {
 tinykepler_state target_at_arrival;
 double transfer_departure_m_s[3],transfer_arrival_m_s[3];
 double departure_delta_v_m_s[3],arrival_delta_v_m_s[3],total_delta_v_m_s;
} tinyrendezvous_intercept;
/* Coast target with two-body Kepler, then solve the selected zero-revolution
   Lambert transfer from chaser's current position to target arrival position.
   Inputs are same-epoch, same inertial frame, SI units. Flight time>0. The
   optional arrival matching impulse is included in total delta-v. */
tinyrendezvous_status tinyrendezvous_target(double mu_m3_s2,const tinykepler_state *chaser,const tinykepler_state *target,double flight_time_s,tinylambert_path path,tinyrendezvous_intercept *out);
typedef struct {double characteristic_delta_v_m_s,ideal_duration_s;} tinyrendezvous_low_thrust;
/* Edelbaum quasi-circular radius/inclination transfer ESTIMATE, constant
   acceleration, no drag/J2/RAAN constraint/eclipses/mass-flow guidance. Bounds:
   |inclination change|<=pi/2 and acceleration<=1% of weaker endpoint gravity. */
tinyrendezvous_status tinyrendezvous_edelbaum(double mu_m3_s2,double initial_radius_m,double final_radius_m,double inclination_change_rad,double acceleration_m_s2,tinyrendezvous_low_thrust *out);
#ifdef __cplusplus
}
#endif
#endif
