#include "tinydisturbance/tinydisturbance.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL %d: %s\n",__LINE__,#c);return 1;}}while(0)
#include "data/reference.h"
static int same(const double a[3],const double b[3]){for(size_t i=0;i<3;++i)if(fabs(a[i]-b[i])>3e-12*fmax(1e-15,fabs(b[i])))return 0;return 1;}
int main(void){double force[3],torque[3],maximum=0;
 for(size_t c=0;c<sizeof reference_cases/sizeof reference_cases[0];++c){CHECK(tinydisturbance_gravity_gradient(3.986004418e14,reference_cases[c].position,reference_cases[c].inertia,torque)==TINYDISTURBANCE_OK);CHECK(same(torque,reference_cases[c].gravity));CHECK(tinydisturbance_dipole(reference_cases[c].dipole,reference_cases[c].field,torque)==TINYDISTURBANCE_OK);CHECK(same(torque,reference_cases[c].magnetic));CHECK(tinydisturbance_aerodynamic(1e-12,2.2,.03,reference_cases[c].normal,reference_cases[c].velocity,reference_cases[c].arm,force,torque)==TINYDISTURBANCE_OK);CHECK(same(force,reference_cases[c].aeroforce)&&same(torque,reference_cases[c].aerotorque));double power=0;for(size_t i=0;i<3;++i)power+=force[i]*reference_cases[c].velocity[i];CHECK(power<=1e-20);CHECK(tinydisturbance_solar_pressure(4.56e-6,.03,.3,.2,reference_cases[c].normal,reference_cases[c].sun,reference_cases[c].arm,force,torque)==TINYDISTURBANCE_OK);CHECK(same(force,reference_cases[c].solarforce)&&same(torque,reference_cases[c].solartorque));for(size_t i=0;i<3;++i)maximum=fmax(maximum,fabs(torque[i]-reference_cases[c].solartorque[i]));}
 {double n[]={1,0,0},sun[]={-1,0,0},arm[]={0,0,1};CHECK(tinydisturbance_solar_pressure(1,1,0,0,n,sun,arm,force,torque)==TINYDISTURBANCE_OK);CHECK(force[0]==0&&torque[1]==0);sun[0]=1;CHECK(tinydisturbance_solar_pressure(1,1,1,0,n,sun,arm,force,torque)==TINYDISTURBANCE_OK);CHECK(force[0]==-2&&torque[1]==-2);}
 {double zero[]={0,0,0},j[]={1,0,0,0,1,0,0,0,1},saved[]={123,456,789};memcpy(torque,saved,sizeof torque);CHECK(tinydisturbance_gravity_gradient(1,zero,j,torque)==TINYDISTURBANCE_DEGENERATE);CHECK(!memcmp(torque,saved,sizeof torque));CHECK(tinydisturbance_solar_pressure(1,1,.8,.5,zero,zero,zero,force,torque)==TINYDISTURBANCE_INVALID_ARGUMENT);}
 printf("tinydisturbance: 64 NumPy multi-model cases; max solar-torque error %.3g N m\n",maximum);return 0;}
