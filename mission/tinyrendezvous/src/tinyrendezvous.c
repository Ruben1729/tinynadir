/* SPDX-License-Identifier: MIT; independent orbital phasing and composition. */
#include "tinyrendezvous/tinyrendezvous.h"
#include <math.h>
#define PI 3.141592653589793238462643383279502884
static double norm(const double v[3]){return hypot(hypot(v[0],v[1]),v[2]);}
tinyrendezvous_status tinyrendezvous_phase_circular(double mu,double radius,double phase,unsigned target_turns,unsigned chaser_turns,double minimum,tinyrendezvous_phasing *out) {
    tinyrendezvous_phasing r;double circular,n,phase_speed;
    if(!out||!isfinite(mu)||!isfinite(radius)||!isfinite(phase)||!isfinite(minimum))return TINYRENDEZVOUS_INVALID;
    if(mu<=0||radius<=0||fabs(phase)>PI||minimum<0||radius<minimum||target_turns<1||chaser_turns<1||target_turns>100000||chaser_turns>100000)return TINYRENDEZVOUS_RANGE;
    circular=sqrt(mu/radius);n=circular/radius;
    if(!isfinite(circular)||n<=0||!isfinite(n))return TINYRENDEZVOUS_RANGE;
    r.duration_s=(2*PI*target_turns-phase)/n;r.phasing_period_s=r.duration_s/chaser_turns;
    r.semimajor_axis_m=cbrt(mu)*pow(r.phasing_period_s/(2*PI),2./3.);
    r.other_apsis_radius_m=2*r.semimajor_axis_m-radius;
    if(!isfinite(r.duration_s)||!isfinite(r.semimajor_axis_m)||r.other_apsis_radius_m<=0||r.other_apsis_radius_m<minimum)return TINYRENDEZVOUS_GEOMETRY;
    phase_speed=sqrt(mu/radius*(2-radius/r.semimajor_axis_m));
    r.departure_tangent_delta_v_m_s=phase_speed-circular;r.arrival_tangent_delta_v_m_s=-r.departure_tangent_delta_v_m_s;
    r.total_delta_v_m_s=2*fabs(r.departure_tangent_delta_v_m_s);
    if(!isfinite(r.total_delta_v_m_s))return TINYRENDEZVOUS_RANGE;
    *out=r;return TINYRENDEZVOUS_OK;
}
tinyrendezvous_status tinyrendezvous_target(double mu,const tinykepler_state *chaser,const tinykepler_state *target,double dt,tinylambert_path path,tinyrendezvous_intercept *out) {
    tinyrendezvous_intercept r;tinylambert_solution solution;int i;
    if(!out||!chaser||!target||!isfinite(mu)||!isfinite(dt))return TINYRENDEZVOUS_INVALID;
    if(mu<=0||dt<=0||(path!=TINYLAMBERT_SHORT_WAY&&path!=TINYLAMBERT_LONG_WAY))return TINYRENDEZVOUS_RANGE;
    for(i=0;i<3;i++)if(!isfinite(chaser->position_m[i])||!isfinite(chaser->velocity_m_s[i])||!isfinite(target->position_m[i])||!isfinite(target->velocity_m_s[i]))return TINYRENDEZVOUS_INVALID;
    if(tinykepler_propagate(mu,target,dt,&r.target_at_arrival)!=TINYKEPLER_OK)return TINYRENDEZVOUS_PROPAGATION;
    if(tinylambert_solve(mu,chaser->position_m,r.target_at_arrival.position_m,dt,path,&solution)!=TINYLAMBERT_OK)return TINYRENDEZVOUS_LAMBERT;
    for(i=0;i<3;i++) {
        r.transfer_departure_m_s[i]=solution.departure_m_s[i];r.transfer_arrival_m_s[i]=solution.arrival_m_s[i];
        r.departure_delta_v_m_s[i]=solution.departure_m_s[i]-chaser->velocity_m_s[i];
        r.arrival_delta_v_m_s[i]=r.target_at_arrival.velocity_m_s[i]-solution.arrival_m_s[i];
    }
    r.total_delta_v_m_s=norm(r.departure_delta_v_m_s)+norm(r.arrival_delta_v_m_s);
    if(!isfinite(r.total_delta_v_m_s))return TINYRENDEZVOUS_RANGE;
    *out=r;return TINYRENDEZVOUS_OK;
}
tinyrendezvous_status tinyrendezvous_edelbaum(double mu,double r0,double r1,double inc,double acceleration,tinyrendezvous_low_thrust *out) {
    tinyrendezvous_low_thrust r;double v0,v1,gravity;
    if(!out||!isfinite(mu)||!isfinite(r0)||!isfinite(r1)||!isfinite(inc)||!isfinite(acceleration))return TINYRENDEZVOUS_INVALID;
    if(mu<=0||r0<=0||r1<=0||fabs(inc)>PI/2||acceleration<=0)return TINYRENDEZVOUS_RANGE;
    v0=sqrt(mu/r0);v1=sqrt(mu/r1);gravity=fmin(mu/r0/r0,mu/r1/r1);
    if(!isfinite(v0)||!isfinite(v1)||gravity<=0||acceleration>.01*gravity)return TINYRENDEZVOUS_RANGE;
    r.characteristic_delta_v_m_s=hypot(v0-v1,2*sqrt(v0)*sqrt(v1)*sin(PI*fabs(inc)/4));
    r.ideal_duration_s=r.characteristic_delta_v_m_s/acceleration;
    if(!isfinite(r.characteristic_delta_v_m_s)||!isfinite(r.ideal_duration_s))return TINYRENDEZVOUS_RANGE;
    *out=r;return TINYRENDEZVOUS_OK;
}
