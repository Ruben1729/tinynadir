/* SPDX-License-Identifier: MIT */
#include "tinypd/tinypd.h"
#include "tinywheels/tinywheels.h"
#include <math.h>
#include <stdio.h>
int main(void) {
    /* An external estimator and guidance supply these two attitudes. */
    const double current[4] = {0.9987502603949663, 0.04997916927067833, 0, 0};
    const double desired[4] = {1,0,0,0}, rate_error[3] = {0,0,0};
    const double kp[3] = {.001,.001,.001}, kd[3] = {.002,.002,.002};
    const double axes[9] = {1,0,0,0,1,0,0,0,1};
    const double speed[3] = {0,0,0}, torque_limit[3] = {.001,.001,.001};
    const double speed_limit[3] = {1000,1000,1000};
    tinywheels_workspace workspace;
    double torque[3], wheel_torque[3], residual[3];
    if (tinypd_torque(current,desired,rate_error,kp,kd,NULL,torque)) return 1;
    if (tinywheels_allocate(3,axes,torque,speed,torque_limit,speed_limit,NULL,wheel_torque,residual,&workspace)) return 1;
    for (size_t i=0;i<3;++i) if (fabs(residual[i])>1e-12) return 1;
    printf("Wheel rotor torques [N m]: %.9g %.9g %.9g\n",wheel_torque[0],wheel_torque[1],wheel_torque[2]);
    return 0;
}
