/* SPDX-License-Identifier: MIT */
#ifndef TINYHILLS_H
#define TINYHILLS_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum { TINYHILLS_OK, TINYHILLS_INVALID, TINYHILLS_SINGULAR } tinyhills_status;
typedef struct { double position_m[3], velocity_m_s[3]; } tinyhills_state;
/* Linear Clohessy-Wiltshire dynamics: x radial outward, y along-track,
   z orbit normal. Chief circular mean motion n in rad/s, elapsed time seconds.
   n>=0, |n*dt|<=1e8. n=0 gives constant-velocity Cartesian motion.
   Assumes separation small relative to chief radius. Failure preserves output. */
tinyhills_status tinyhills_propagate(double mean_motion_rad_s,double elapsed_s,
                                   const tinyhills_state *initial,tinyhills_state *output);
/* Row-major 6x6 state transition, state ordered [x,y,z,vx,vy,vz]. */
tinyhills_status tinyhills_transition(double mean_motion_rad_s,double elapsed_s,double transition[36]);
/* Initial relative velocity to arrive at relative position zero after dt.
   Returns SINGULAR when the terminal-position map is not uniquely invertible. */
tinyhills_status tinyhills_rendezvous(double mean_motion_rad_s,double elapsed_s,
                                    const double initial_position_m[3],double velocity_m_s[3]);
#ifdef __cplusplus
}
#endif
#endif
