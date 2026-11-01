/* SPDX-License-Identifier: MIT */
#ifndef TINYKEPLER_H
#define TINYKEPLER_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum { TINYKEPLER_OK, TINYKEPLER_INVALID, TINYKEPLER_RANGE, TINYKEPLER_NO_CONVERGENCE } tinykepler_status;
typedef struct { double position_m[3], velocity_m_s[3]; } tinykepler_state;
/* Osculating nonradial state in any one inertial frame, mu in m^3/s^2.
   Elliptic, parabolic, hyperbolic two-body motion; no perturbations.
   |nondimensional elapsed time| <= 1e8 after elliptic period reduction.
   Bracketing and root solving are bounded. Input/output may alias.
   Failure preserves output. */
tinykepler_status tinykepler_propagate(double mu_m3_s2,
                                     const tinykepler_state *initial,
                                     double elapsed_s, tinykepler_state *output);
#ifdef __cplusplus
}
#endif
#endif
