/* SPDX-License-Identifier: MIT */
#ifndef TINYHOHMANN_H
#define TINYHOHMANN_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum { TINYHOHMANN_OK, TINYHOHMANN_INVALID, TINYHOHMANN_GEOMETRY } tinyhohmann_status;
typedef struct { double burn_m_s[3],total_delta_v_m_s,time_of_flight_s;unsigned burn_count; } tinyhohmann_transfer;
typedef struct { tinyhohmann_transfer transfer;double eccentricity,p_m,arrival_flight_path_rad; } tinyhohmann_tangent;
/* Coplanar impulsive transfers between circular orbits. Radii metres, mu m^3/s^2.
   Burns are nonnegative magnitudes, not inertial vectors or actuator commands.
   Bielliptic intermediate radius must be >= both endpoint radii.
   Failures preserve output. */
tinyhohmann_status tinyhohmann_circular(double mu,double initial_radius_m,double final_radius_m,tinyhohmann_transfer *output);
tinyhohmann_status tinyhohmann_bielliptic(double mu,double initial_radius_m,double intermediate_radius_m,double final_radius_m,tinyhohmann_transfer *output);
/* One-tangent coast between circular endpoint orbits. Starts at transfer
   periapsis (start_at_apoapsis=0, final true anomaly in (0,pi]) or apoapsis
   (nonzero, final true anomaly in (pi,2*pi)). The physical geometry determines
   transfer eccentricity; impossible/nonphysical transfers are rejected. */
tinyhohmann_status tinyhohmann_one_tangent(double mu,double initial_radius_m,double final_radius_m,
                                         int start_at_apoapsis,double final_true_anomaly_rad,tinyhohmann_tangent *output);
#ifdef __cplusplus
}
#endif
#endif
