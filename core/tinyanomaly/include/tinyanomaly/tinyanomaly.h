/* SPDX-License-Identifier: MIT */
#ifndef TINYANOMALY_H
#define TINYANOMALY_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum { TINYANOMALY_OK, TINYANOMALY_INVALID, TINYANOMALY_NO_CONVERGENCE } tinyanomaly_status;
typedef struct { double anomaly_rad; double true_anomaly_rad; } tinyanomaly_solution;
/* e in [0,1000], |mean|<=1e12. Elliptic mean is reduced to [-pi,pi].
   anomaly_rad is eccentric E for e<1, hyperbolic H for e>1, Barker D for e=1.
   For e=1 the mean argument is D+D^3/3, dimensionless (not an angle).
   Both functions preserve outputs on failure. All angles in radians. */
tinyanomaly_status tinyanomaly_from_mean(double eccentricity, double mean,
                                  tinyanomaly_solution *output);
tinyanomaly_status tinyanomaly_from_true(double eccentricity, double true_anomaly_rad,
                                  double *mean, double *anomaly);
#ifdef __cplusplus
}
#endif
#endif
