#include "tinyrendezvous/tinyrendezvous.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(x) do{if(!(x)){fprintf(stderr,"%s:%d %s\n",__FILE__,__LINE__,#x);exit(1);}}while(0)
#define PI 3.14159265358979323846
#include "data/reference.h"
int main(void) {
 size_t i;int k;tinyrendezvous_phasing phase,saved;tinyrendezvous_low_thrust low;tinyrendezvous_intercept intercept,before;
 tinykepler_state chaser,target,coast;double maxphase=0,maxlow=0,maxv=0;
 for(i=0;i<sizeof phasing_cases/sizeof phasing_cases[0];i++) {
  const double *q=phasing_cases[i];double err,n;
  CHECK(tinyrendezvous_phase_circular(q[0],q[1],q[2],(unsigned)q[3],(unsigned)q[4],q[5],&phase)==TINYRENDEZVOUS_OK);
  CHECK(fabs(phase.duration_s-q[6])<1e-7);err=fabs(phase.total_delta_v_m_s-q[7]);CHECK(err<1e-7);if(err>maxphase)maxphase=err;
  n=sqrt(q[0]/q[1])/q[1];CHECK(fabs(remainder(q[2]+n*phase.duration_s,2*PI))<1e-13);
  CHECK(fabs(phase.phasing_period_s*q[4]-phase.duration_s)<1e-9);
  CHECK(fabs(phase.departure_tangent_delta_v_m_s+phase.arrival_tangent_delta_v_m_s)<1e-12);
 }
 for(i=0;i<sizeof edelbaum_cases/sizeof edelbaum_cases[0];i++) {
  const double *q=edelbaum_cases[i];double error;
  CHECK(tinyrendezvous_edelbaum(q[0],q[1],q[2],q[3],q[4],&low)==TINYRENDEZVOUS_OK);
  error=fabs(low.characteristic_delta_v_m_s-q[5]);CHECK(error<1e-7);if(error>maxlow)maxlow=error;CHECK(fabs(low.ideal_duration_s-q[6])<.001);
 }
 for(i=0;i<sizeof intercept_cases/sizeof intercept_cases[0];i++) {
  const double *q=intercept_cases[i];double error;
  for(k=0;k<3;k++){chaser.position_m[k]=q[1+k];chaser.velocity_m_s[k]=q[4+k];target.position_m[k]=q[7+k];target.velocity_m_s[k]=q[10+k];}
  CHECK(tinyrendezvous_target(q[0],&chaser,&target,q[13],(tinylambert_path)(int)q[14],&intercept)==TINYRENDEZVOUS_OK);
  for(k=0;k<3;k++){
   CHECK(fabs(intercept.target_at_arrival.position_m[k]-q[15+k])<.001);CHECK(fabs(intercept.target_at_arrival.velocity_m_s[k]-q[18+k])<1e-6);
   error=fabs(intercept.transfer_departure_m_s[k]-q[21+k]);CHECK(error<1e-6);if(error>maxv)maxv=error;
   CHECK(fabs(intercept.transfer_arrival_m_s[k]-q[24+k])<1e-6);CHECK(fabs(intercept.departure_delta_v_m_s[k]-q[27+k])<1e-6);CHECK(fabs(intercept.arrival_delta_v_m_s[k]-q[30+k])<1e-6);
  }
  CHECK(fabs(intercept.total_delta_v_m_s-q[33])<1e-5);
  coast=chaser;for(k=0;k<3;k++)coast.velocity_m_s[k]=intercept.transfer_departure_m_s[k];
  CHECK(tinykepler_propagate(q[0],&coast,q[13],&coast)==TINYKEPLER_OK);
  for(k=0;k<3;k++)CHECK(fabs(coast.position_m[k]-intercept.target_at_arrival.position_m[k])<.001);
 }
 saved=phase;CHECK(tinyrendezvous_phase_circular(3.986e14,7e6,1,1,1,6.9e6,&phase)==TINYRENDEZVOUS_GEOMETRY);CHECK(!memcmp(&phase,&saved,sizeof phase));
 CHECK(tinyrendezvous_phase_circular(3.986e14,7e6,NAN,1,1,0,&phase)==TINYRENDEZVOUS_INVALID);
 CHECK(tinyrendezvous_edelbaum(3.986e14,7e6,8e6,.1,20,&low)==TINYRENDEZVOUS_RANGE);
 before=intercept;CHECK(tinyrendezvous_target(3.986e14,&chaser,&target,0,TINYLAMBERT_SHORT_WAY,&intercept)==TINYRENDEZVOUS_RANGE);CHECK(!memcmp(&before,&intercept,sizeof intercept));
 printf("tinyrendezvous: 45 Vallado phasing, 20 poliastro Edelbaum, %zu moving-target cases pass; maxima %.9g / %.9g / %.9g m/s\n",sizeof intercept_cases/sizeof intercept_cases[0],maxphase,maxlow,maxv);return 0;
}
