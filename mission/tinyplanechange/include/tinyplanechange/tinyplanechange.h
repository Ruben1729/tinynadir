/* SPDX-License-Identifier: MIT */
#ifndef TINYPLANECHANGE_H
#define TINYPLANECHANGE_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum { TINYPLANECHANGE_OK,TINYPLANECHANGE_INVALID } tinyplanechange_status;
typedef struct {double delta_v_m_s,plane_angle_rad,initial_argument_latitude_rad,final_argument_latitude_rad;int unique_intersection;} tinyplanechange_result;
/* Ideal instantaneous same-speed plane change. Angles radians, speed m/s.
   The radial velocity component is unchanged; only v*cos(flight_path_angle)
   turns. Inclinations in [0,pi]; flight path in [-pi/2,pi/2]. Failures preserve outputs. */
tinyplanechange_status tinyplanechange_inclination(double speed_m_s,double flight_path_rad,double delta_inclination_rad,double *delta_v_m_s);
tinyplanechange_status tinyplanechange_combined(double speed_m_s,double flight_path_rad,double initial_inclination_rad,double final_inclination_rad,double delta_raan_rad,tinyplanechange_result *output);
tinyplanechange_status tinyplanechange_node(double speed_m_s,double flight_path_rad,double inclination_rad,double delta_raan_rad,tinyplanechange_result *output);
#ifdef __cplusplus
}
#endif
#endif
