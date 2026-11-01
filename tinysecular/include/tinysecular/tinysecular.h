/* SPDX-License-Identifier: MIT */
#ifndef TINYSECULAR_H
#define TINYSECULAR_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYSECULAR_OK=0,TINYSECULAR_INVALID,TINYSECULAR_RANGE} tinysecular_status;
typedef struct {double mu_m3_s2,equatorial_radius_m,j2;} tinysecular_model;
typedef struct {double semimajor_axis_m,eccentricity,inclination_rad,raan_rad,arg_periapsis_rad,mean_anomaly_rad;} tinysecular_elements;
/* Both profiles use the same first-order J2 nodal/apsidal precession.
   AVERAGED_J2 includes the first-order mean-anomaly secular correction.
   VALLADO_PKEPLER holds mean motion at sqrt(mu/a^3), as Vallado Algorithm66.
   These are explicitly defined mean elements, never SGP4/TLE mean elements. */
typedef enum {TINYSECULAR_AVERAGED_J2=0,TINYSECULAR_VALLADO_PKEPLER} tinysecular_profile;
typedef struct {double raan_rad_s,arg_periapsis_rad_s,mean_anomaly_rad_s;} tinysecular_rates;
/* Elliptic orbits with pericenter strictly outside the model radius; |J2|<=.01.
   Angles radians, SI otherwise; axial gravity pole is the inertial +z axis.
   No short-period terms, drag, J3/J4 or third bodies. Outputs atomic on error. */
tinysecular_status tinysecular_j2_rates(const tinysecular_model *model,const tinysecular_elements *mean,tinysecular_profile profile,tinysecular_rates *out);
/* a,e,i held fixed; angles wrapped to [0,2pi). |rate*time| limited to 1e12 rad. */
tinysecular_status tinysecular_propagate(const tinysecular_model *model,const tinysecular_elements *initial,tinysecular_profile profile,double elapsed_s,tinysecular_elements *out);
#ifdef __cplusplus
}
#endif
#endif
