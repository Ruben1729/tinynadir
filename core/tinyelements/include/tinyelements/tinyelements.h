/* SPDX-License-Identifier: MIT */
#ifndef TINYELEMENTS_H
#define TINYELEMENTS_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum { TINYELEMENTS_OK, TINYELEMENTS_INVALID, TINYELEMENTS_SINGULAR } tinyelements_status;
typedef struct { double position_m[3], velocity_m_s[3]; } tinyelements_state;
typedef struct {
    double p_m,eccentricity,inclination_rad,raan_rad,arg_periapsis_rad,true_anomaly_rad;
} tinyelements_classical;
/* Modified equinoctial elements: f=e*cos(Omega+omega), g=e*sin(Omega+omega),
   h=tan(i/2)*cos(Omega), k=tan(i/2)*sin(Omega), L=Omega+omega+nu.
   Prograde convention, singular at i=pi. */
typedef struct { double p_m,f,g,h,k,longitude_rad; } tinyelements_equinoctial;
typedef enum { TINYELEMENTS_RSW, TINYELEMENTS_NTW, TINYELEMENTS_PQW } tinyelements_frame;
/* SI units and radians, mu in m^3/s^2. Angles returned in [0,2*pi).
   e<1e-12: e=0, omega=0 and nu=argument of latitude.
   |sin(i)|<1e-12: Omega=0, omega is signed longitude of periapsis.
   Circular equatorial: Omega=omega=0, nu is signed true longitude.
   Nonradial elliptic/parabolic/hyperbolic states supported; failures preserve output. */
tinyelements_status tinyelements_from_classical(double mu,const tinyelements_classical *elements,tinyelements_state *state);
tinyelements_status tinyelements_to_classical(double mu,const tinyelements_state *state,tinyelements_classical *elements);
tinyelements_status tinyelements_from_equinoctial(double mu,const tinyelements_equinoctial *elements,tinyelements_state *state);
tinyelements_status tinyelements_to_equinoctial(double mu,const tinyelements_state *state,tinyelements_equinoctial *elements);
/* Row-major basis maps inertial vector coordinates to RSW/NTW/PQW coordinates.
   This is a vector rotation, not a rotating-frame relative-velocity conversion. */
tinyelements_status tinyelements_local_frame(double mu,const tinyelements_state *state,tinyelements_frame frame,double rotation[9]);
#ifdef __cplusplus
}
#endif
#endif
