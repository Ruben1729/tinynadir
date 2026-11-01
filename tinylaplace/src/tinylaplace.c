/* SPDX-License-Identifier: MIT */
#include "tinylaplace/tinylaplace.h"
#include "tinylinalg/tinylinalg.h"
#include <math.h>
#include <stddef.h>
#include <string.h>
static int finite_array(const double *v,size_t n){if(!v)return 0;for(size_t i=0;i<n;++i)if(!isfinite(v[i]))return 0;return 1;}
static double triple(const double a[3],const double b[3],const double c[3]){return a[0]*(b[1]*c[2]-b[2]*c[1])+a[1]*(b[2]*c[0]-b[0]*c[2])+a[2]*(b[0]*c[1]-b[1]*c[0]);}
tinylaplace_status tinylaplace_estimate(double mu,const double l[3],const double ld[3],const double ldd[3],const double position[3],
    const double velocity[3],const double acceleration[3],double hint,tinylaplace_result *result,tinylaplace_workspace *work){
    double scale,R[3],V[3],acc[3],d,A,B,projection,mu_scaled,coefficients[9]={1,0,0,0,0,0,0,0,0},best_score=HUGE_VAL;
    tinylaplace_result best={0};unsigned candidates=0;
    if(!result||!work||!isfinite(mu)||mu<=0||!isfinite(hint)||hint<0||!finite_array(l,3)||!finite_array(ld,3)||!finite_array(ldd,3)||!finite_array(position,3)||!finite_array(velocity,3)||!finite_array(acceleration,3))return TINYLAPLACE_INVALID_ARGUMENT;
    if(fabs(tinylinalg_norm3(l)-1)>1e-10)return TINYLAPLACE_INVALID_ARGUMENT;
    scale=tinylinalg_norm3(position);if(!isfinite(scale)||scale==0)return TINYLAPLACE_GEOMETRY;
    d=triple(l,ld,ldd);if(!isfinite(d)||fabs(d)<=1e-12*tinylinalg_norm3(ld)*tinylinalg_norm3(ldd))return TINYLAPLACE_GEOMETRY;
    for(size_t i=0;i<3;++i){R[i]=position[i]/scale;V[i]=velocity[i]/scale;acc[i]=acceleration[i]/scale;}
    mu_scaled=((mu/scale)/scale)/scale;A=-triple(l,ld,acc)/d;B=-mu_scaled*triple(l,ld,R)/d;projection=tinylinalg_dot3(l,R);
    coefficients[2]=-(1+A*A+2*A*projection);coefficients[5]=-2*B*(A+projection);coefficients[8]=-B*B;
    if(!finite_array(coefficients,9))return TINYLAPLACE_ARITHMETIC_ERROR;
    if(tinypolynomial_roots(8,coefficients,work->roots,work->scratch,16)!=TINYPOLYNOMIAL_OK)return TINYLAPLACE_NO_CONVERGENCE;
    for(size_t root=0;root<8;++root){
        double radius=work->roots[root].real,u,rho,rhodot,score;tinylaplace_result candidate={0};
        if(radius<=0||fabs(work->roots[root].imaginary)>1e-8*fmax(1,fabs(radius)))continue;
        u=mu_scaled/radius/radius/radius;rho=A+B/radius/radius/radius;
        if(rho<=0||!isfinite(rho))continue;
        rhodot=(triple(l,ldd,acc)+u*triple(l,ldd,R))/(2*d);
        candidate.range_m=rho*scale;candidate.range_rate_m_s=rhodot*scale;candidate.selected_radius_m=radius*scale;
        for(size_t i=0;i<3;++i){candidate.position_m[i]=(R[i]+rho*l[i])*scale;candidate.velocity_m_s[i]=(V[i]+rhodot*l[i]+rho*ld[i])*scale;}
        if(!finite_array(candidate.position_m,3)||!finite_array(candidate.velocity_m_s,3)||!isfinite(candidate.selected_radius_m)||!isfinite(candidate.range_m)||!isfinite(candidate.range_rate_m_s))continue;
        ++candidates;score=hint>0?fabs(radius-hint/scale):-radius;if(score<best_score){best=candidate;best_score=score;}
    }
    if(!candidates)return TINYLAPLACE_NO_SOLUTION;
    best.admissible_roots=candidates;*result=best;return TINYLAPLACE_OK;
}
tinylaplace_status tinylaplace_three(double mu,double spacing,const double los[9],const double position[3],const double velocity[3],
    const double acceleration[3],double hint,tinylaplace_result *result,tinylaplace_workspace *work){
    double l[9],ld[3],ldd[3];
    if(!isfinite(spacing)||spacing<=0||spacing>1e6||!finite_array(los,9))return TINYLAPLACE_INVALID_ARGUMENT;
    for(size_t i=0;i<3;++i)if(tinylinalg_unit3(los+3*i,l+3*i)!=TINYLINALG_OK)return TINYLAPLACE_GEOMETRY;
    for(size_t i=0;i<3;++i){ld[i]=(l[6+i]-l[i])/(2*spacing);ldd[i]=((l[6+i]-l[3+i])+(l[i]-l[3+i]))/spacing/spacing;}
    return tinylaplace_estimate(mu,l+3,ld,ldd,position,velocity,acceleration,hint,result,work);
}
