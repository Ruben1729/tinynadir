/* SPDX-License-Identifier: MIT */
#include "tinykepler/tinykepler.h"
#include "tinystumpff/tinystumpff.h"
#include <math.h>
#include <float.h>

static int equation(double x,double beta,double radial,double tau,
                    double *value,double *derivative,tinystumpff_values *c)
{
    double x2=x*x;
    if (!isfinite(x2) || tinystumpff_evaluate(beta*x2,c)!=TINYSTUMPFF_OK) return 0;
    *value=radial*x2*c->c[2]+(1.0-beta)*x2*x*c->c[3]+x-tau;
    *derivative=1.0+radial*x*c->c[1]+(1.0-beta)*x2*c->c[2];
    return isfinite(*value)&&isfinite(*derivative)&&*derivative>0.0;
}

tinykepler_status tinykepler_propagate(double mu,const tinykepler_state *initial,
                                     double dt,tinykepler_state *output)
{
    double r0,scale_v,scale_t,r[3],v[3],v2=0.0,radial=0.0,h[3],hmag;
    double beta,tau,lo,hi,x,F,D,period,f,g,rmag,fdot,gdot;
    tinystumpff_values c;
    tinykepler_state result;
    unsigned i,iteration;
    if (!initial || !output || !isfinite(mu) || mu<=0.0 || !isfinite(dt)) return TINYKEPLER_INVALID;
    for(i=0;i<3;++i)
        if (!isfinite(initial->position_m[i]) || !isfinite(initial->velocity_m_s[i])) return TINYKEPLER_INVALID;
    r0=hypot(hypot(initial->position_m[0],initial->position_m[1]),initial->position_m[2]);
    if (!(r0>0.0) || !isfinite(r0)) return TINYKEPLER_INVALID;
    scale_v=sqrt(mu/r0);scale_t=r0/scale_v;
    if (!(scale_v>0.0) || !isfinite(scale_v) || !isfinite(scale_t)) return TINYKEPLER_RANGE;
    for(i=0;i<3;++i) {
        r[i]=initial->position_m[i]/r0;v[i]=initial->velocity_m_s[i]/scale_v;
        radial+=r[i]*v[i];v2+=v[i]*v[i];
    }
    h[0]=r[1]*v[2]-r[2]*v[1];h[1]=r[2]*v[0]-r[0]*v[2];h[2]=r[0]*v[1]-r[1]*v[0];
    hmag=hypot(hypot(h[0],h[1]),h[2]);
    if (!isfinite(v2) || !isfinite(hmag)) return TINYKEPLER_RANGE;
    if (hmag<=1e-12*fmax(1.0,sqrt(v2))) return TINYKEPLER_INVALID;
    beta=2.0-v2;tau=dt/scale_t;
    if (!isfinite(tau)) return TINYKEPLER_RANGE;
    if (beta>0.0) {
        period=6.283185307179586476925286766559/(beta*sqrt(beta));
        if (isfinite(period)) tau=remainder(tau,period);
    }
    if (fabs(tau)>1e8) return TINYKEPLER_RANGE;
    if (tau==0.0) { *output=*initial;return TINYKEPLER_OK; }
    lo=tau>0.0 ? 0.0 : -1.0;hi=tau>0.0 ? 1.0 : 0.0;
    for(iteration=0;iteration<128;++iteration) {
        double bound=tau>0.0 ? hi : lo;
        if (!equation(bound,beta,radial,tau,&F,&D,&c)) return TINYKEPLER_RANGE;
        if ((tau>0.0 && F>=0.0) || (tau<0.0 && F<=0.0)) break;
        if (tau>0.0) hi*=2.0; else lo*=2.0;
    }
    if (iteration==128) return TINYKEPLER_NO_CONVERGENCE;
    x=0.5*(lo+hi);
    for(iteration=0;iteration<128;++iteration) {
        double next;
        if (!equation(x,beta,radial,tau,&F,&D,&c)) return TINYKEPLER_RANGE;
        if (F==0.0) break;
        if (F>0.0) hi=x; else lo=x;
        next=x-F/D;
        if (!(next>lo && next<hi)) next=lo+0.5*(hi-lo);
        if (next==x || hi-lo<=4.0*DBL_EPSILON*(fabs(x)+1e-8)) break;
        x=next;
    }
    if (iteration==128) return TINYKEPLER_NO_CONVERGENCE;
    if (!equation(x,beta,radial,tau,&F,&D,&c)) return TINYKEPLER_RANGE;
    f=1.0-x*x*c.c[2];g=tau-x*x*x*c.c[3];
    for(i=0;i<3;++i) result.position_m[i]=f*r[i]+g*v[i];
    rmag=hypot(hypot(result.position_m[0],result.position_m[1]),result.position_m[2]);
    if (!(rmag>0.0) || !isfinite(rmag)) return TINYKEPLER_RANGE;
    fdot=(beta*x*x*x*c.c[3]-x)/rmag;gdot=1.0-x*x*c.c[2]/rmag;
    for(i=0;i<3;++i) {
        result.velocity_m_s[i]=(fdot*r[i]+gdot*v[i])*scale_v;
        result.position_m[i]*=r0;
        if (!isfinite(result.position_m[i]) || !isfinite(result.velocity_m_s[i])) return TINYKEPLER_RANGE;
    }
    *output=result;
    return TINYKEPLER_OK;
}
