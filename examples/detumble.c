/* SPDX-License-Identifier: MIT */
#include "tinybdot/tinybdot.h"
#include <stdio.h>
int main(void) {
    tinybdot_state state = {0};
    const double samples[2][3] = {{20e-6, -5e-6, 35e-6}, {19e-6, -4e-6, 34e-6}};
    const double limit[3] = {.02, .02, .02};
    double dipole[3];
    unsigned limited;
    /* Calibrated body-frame magnetometer data; acquire with torquers quiet.
       No orbit, Sun model, IGRF, attitude solver, or gyro is required here. */
    for (size_t i=0; i<2; ++i) {
        if (tinybdot_step(&state, samples[i], .1, 1000, 1, limit, dipole, &limited)) return 1;
    }
    printf("B-dot dipole [A m^2]: %.9g %.9g %.9g; limited axes %u\n", dipole[0],dipole[1],dipole[2],limited);
    return 0;
}
