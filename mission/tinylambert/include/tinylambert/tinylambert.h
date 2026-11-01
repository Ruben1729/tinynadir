/* SPDX-License-Identifier: MIT */
#ifndef TINYLAMBERT_H
#define TINYLAMBERT_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum { TINYLAMBERT_OK,TINYLAMBERT_INVALID,TINYLAMBERT_GEOMETRY,TINYLAMBERT_RANGE,TINYLAMBERT_NO_CONVERGENCE } tinylambert_status;
typedef enum { TINYLAMBERT_SHORT_WAY,TINYLAMBERT_LONG_WAY } tinylambert_path;
typedef struct { double departure_m_s[3],arrival_m_s[3];unsigned iterations; } tinylambert_solution;
/* Zero-revolution universal-variable Lambert solution. SI units, inertial
   positions in the same frame, mu m^3/s^2, strictly positive flight time.
   SHORT_WAY chooses angle<pi; LONG_WAY chooses angle>pi. The position pair
   determines their transfer-plane direction; no separate inclination/prograde
   flag is implied. Collinear geometry is rejected. No multirevolution branches.
   Nondimensional flight time must be in [1e-6,1e8]. Failure preserves output. */
tinylambert_status tinylambert_solve(double mu,const double departure_m[3],
 const double arrival_m[3],double flight_time_s,tinylambert_path path,tinylambert_solution *output);
#ifdef __cplusplus
}
#endif
#endif
