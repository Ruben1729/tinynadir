#include "tinysecular/tinysecular.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(x) do{if(!(x)){fprintf(stderr,"%s:%d %s\n",__FILE__,__LINE__,#x);exit(1);}}while(0)
#define PI 3.14159265358979323846
#include "data/reference.h"
int main(void){
 size_t i;int k;tinysecular_model m;tinysecular_elements e,out,saved;tinysecular_rates r,before;double maxrate=0,maxangle=0;
 for(i=0;i<sizeof rates_cases/sizeof rates_cases[0];i++){
  const double*q=rates_cases[i];double errors[3];m.mu_m3_s2=q[0];m.equatorial_radius_m=q[1];m.j2=q[2];e.semimajor_axis_m=q[3];e.eccentricity=q[4];e.inclination_rad=q[5];e.raan_rad=q[6];e.arg_periapsis_rad=q[7];e.mean_anomaly_rad=q[8];
  CHECK(tinysecular_j2_rates(&m,&e,TINYSECULAR_AVERAGED_J2,&r)==TINYSECULAR_OK);
  errors[0]=fabs(r.raan_rad_s-q[9]);errors[1]=fabs(r.arg_periapsis_rad_s-q[10]);errors[2]=fabs(r.mean_anomaly_rad_s-q[11]);
  for(k=0;k<3;k++){CHECK(errors[k]<5e-16);if(errors[k]>maxrate)maxrate=errors[k];}
 }
 for(i=0;i<sizeof prop_cases/sizeof prop_cases[0];i++){
  const double*q=prop_cases[i];double errors[3];m.mu_m3_s2=q[0];m.equatorial_radius_m=q[1];m.j2=q[2];e.semimajor_axis_m=q[3];e.eccentricity=q[4];e.inclination_rad=q[5];e.raan_rad=q[6];e.arg_periapsis_rad=q[7];e.mean_anomaly_rad=q[8];
  CHECK(tinysecular_propagate(&m,&e,TINYSECULAR_VALLADO_PKEPLER,q[9],&out)==TINYSECULAR_OK);
  errors[0]=fabs(remainder(out.raan_rad-q[10],2*PI));errors[1]=fabs(remainder(out.arg_periapsis_rad-q[11],2*PI));errors[2]=fabs(remainder(out.mean_anomaly_rad-q[12],2*PI));
  for(k=0;k<3;k++){CHECK(errors[k]<2e-11);if(errors[k]>maxangle)maxangle=errors[k];}
  CHECK(out.semimajor_axis_m==e.semimajor_axis_m&&out.eccentricity==e.eccentricity&&out.inclination_rad==e.inclination_rad);
  CHECK(tinysecular_propagate(&m,&out,TINYSECULAR_VALLADO_PKEPLER,-q[9],&out)==TINYSECULAR_OK);
  CHECK(fabs(remainder(out.mean_anomaly_rad-e.mean_anomaly_rad,2*PI))<1e-12);
 }
 m.j2=0;CHECK(tinysecular_j2_rates(&m,&e,TINYSECULAR_AVERAGED_J2,&r)==TINYSECULAR_OK);CHECK(r.raan_rad_s==0&&r.arg_periapsis_rad_s==0);CHECK(fabs(r.mean_anomaly_rad_s-sqrt(m.mu_m3_s2/e.semimajor_axis_m)/e.semimajor_axis_m)<1e-18);
 e.eccentricity=0;e.inclination_rad=0;CHECK(tinysecular_j2_rates(&m,&e,TINYSECULAR_AVERAGED_J2,&r)==TINYSECULAR_OK);
 before=r;m.j2=.02;CHECK(tinysecular_j2_rates(&m,&e,TINYSECULAR_AVERAGED_J2,&r)==TINYSECULAR_RANGE);CHECK(!memcmp(&r,&before,sizeof r));m.j2=.001;
 e.semimajor_axis_m=m.equatorial_radius_m;CHECK(tinysecular_j2_rates(&m,&e,TINYSECULAR_AVERAGED_J2,&r)==TINYSECULAR_RANGE);
 saved=out;CHECK(tinysecular_propagate(&m,&e,TINYSECULAR_AVERAGED_J2,NAN,&out)==TINYSECULAR_INVALID);CHECK(!memcmp(&out,&saved,sizeof out));
 printf("tinysecular: 24 SciPy Gauss averages and 72 Vallado pkepler cases; max %.9g rad/s, %.9g rad\n",maxrate,maxangle);return 0;
}
