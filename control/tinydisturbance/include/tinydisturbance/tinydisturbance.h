#ifndef TINYDISTURBANCE_H
#define TINYDISTURBANCE_H
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYDISTURBANCE_OK=0,TINYDISTURBANCE_INVALID_ARGUMENT,TINYDISTURBANCE_DEGENERATE,TINYDISTURBANCE_ARITHMETIC_ERROR} tinydisturbance_status;
/* All vectors are in the same body axes; SI units. Outputs preserved on error.
 * Distinct output buffers are required when both force and torque are returned. */
tinydisturbance_status tinydisturbance_gravity_gradient(double mu,const double center_to_spacecraft[3],const double inertia[9],double torque[3]);
tinydisturbance_status tinydisturbance_dipole(const double residual_dipole[3],const double field[3],double torque[3]);
tinydisturbance_status tinydisturbance_force(const double center_of_pressure_arm[3],const double force[3],double torque[3]);
/* One-sided projected-area drag plate; normal points out of the exposed face.
 * Relative velocity is spacecraft relative to atmosphere. Force opposes it. */
tinydisturbance_status tinydisturbance_aerodynamic(double density,double cd,double area,const double normal[3],
    const double relative_velocity[3],const double arm[3],double force[3],double torque[3]);
/* One-sided optical plate, Sun direction points spacecraft→Sun. specular and
 * diffuse reflectivities >=0 and sum <=1. pressure is incident N/m². */
tinydisturbance_status tinydisturbance_solar_pressure(double pressure,double area,double specular,double diffuse,
    const double normal[3],const double sun_direction[3],const double arm[3],double force[3],double torque[3]);
#ifdef __cplusplus
}
#endif
#endif
