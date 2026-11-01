/* SPDX-License-Identifier: MIT */
#include "tinyelements/tinyelements.h"
#include <math.h>
#include <float.h>
#define TE_PI 3.1415926535897932384626433832795
static double norm(const double a[3]) {return hypot(hypot(a[0],a[1]),a[2]);}
static double dot(const double a[3],const double b[3]) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
static void cross(const double a[3],const double b[3],double c[3])
{c[0]=a[1]*b[2]-a[2]*b[1];c[1]=a[2]*b[0]-a[0]*b[2];c[2]=a[0]*b[1]-a[1]*b[0];}
static double wrap(double a) {a=fmod(a,2*TE_PI);return a<0 ? a+2*TE_PI : a;}
static double angle(const double a[3],const double b[3],const double axis[3])
{
    double an=norm(a),bn=norm(b),u[3],v[3],w[3];unsigned j;
    for(j=0;j<3;++j){u[j]=a[j]/an;v[j]=b[j]/bn;}
    cross(u,v,w);return wrap(atan2(dot(w,axis),dot(u,v)));
}
static int finite_elements(const tinyelements_classical *e)
{
    return e && isfinite(e->p_m) && e->p_m>0 && isfinite(e->eccentricity) && e->eccentricity>=0 &&
        isfinite(e->inclination_rad) && e->inclination_rad>=0 && e->inclination_rad<=TE_PI &&
        isfinite(e->raan_rad) && isfinite(e->arg_periapsis_rad) && isfinite(e->true_anomaly_rad);
}
static void basis(const tinyelements_classical *e,double P[3],double Q[3],double W[3])
{
    double O=e->raan_rad,w=e->arg_periapsis_rad,i=e->inclination_rad;
    double cO=cos(O),sO=sin(O),cw=cos(w),sw=sin(w),ci=cos(i),si=sin(i);
    P[0]=cO*cw-sO*sw*ci;P[1]=sO*cw+cO*sw*ci;P[2]=sw*si;
    Q[0]=-cO*sw-sO*cw*ci;Q[1]=-sO*sw+cO*cw*ci;Q[2]=cw*si;
    W[0]=sO*si;W[1]=-cO*si;W[2]=ci;
}
tinyelements_status tinyelements_from_classical(double mu,const tinyelements_classical *e,tinyelements_state *out)
{
    double P[3],Q[3],W[3],c,s,den,r,v;tinyelements_state result;unsigned j;
    if(!out || !isfinite(mu) || mu<=0 || !finite_elements(e)) return TINYELEMENTS_INVALID;
    c=cos(e->true_anomaly_rad);s=sin(e->true_anomaly_rad);den=1+e->eccentricity*c;
    if(!(den>0)) return TINYELEMENTS_INVALID;
    r=e->p_m/den;v=sqrt(mu/e->p_m);basis(e,P,Q,W);
    for(j=0;j<3;++j) {
        result.position_m[j]=r*(c*P[j]+s*Q[j]);
        result.velocity_m_s[j]=v*(-s*P[j]+(e->eccentricity+c)*Q[j]);
        if(!isfinite(result.position_m[j]) || !isfinite(result.velocity_m_s[j])) return TINYELEMENTS_INVALID;
    }
    *out=result;return TINYELEMENTS_OK;
}
tinyelements_status tinyelements_to_classical(double mu,const tinyelements_state *s,tinyelements_classical *out)
{
    double rn,vn,h[3],hn,W[3],N[3],en[3],scaled_h[3],nlen,ecc;
    tinyelements_classical e;unsigned j;
    if(!out || !s || !isfinite(mu) || mu<=0) return TINYELEMENTS_INVALID;
    for(j=0;j<3;++j) if(!isfinite(s->position_m[j]) || !isfinite(s->velocity_m_s[j])) return TINYELEMENTS_INVALID;
    rn=norm(s->position_m);vn=norm(s->velocity_m_s);
    if(!(rn>0) || !isfinite(rn) || !isfinite(vn)) return TINYELEMENTS_INVALID;
    cross(s->position_m,s->velocity_m_s,h);hn=norm(h);
    if(!isfinite(hn) || hn<=64*DBL_EPSILON*rn*vn) return TINYELEMENTS_SINGULAR;
    for(j=0;j<3;++j) {W[j]=h[j]/hn;scaled_h[j]=h[j]/mu;}
    cross(s->velocity_m_s,scaled_h,en);
    for(j=0;j<3;++j) en[j]-=s->position_m[j]/rn;
    ecc=norm(en);nlen=hypot(W[0],W[1]);
    e.p_m=(hn/mu)*hn;e.eccentricity=ecc<1e-12 ? 0 : ecc;
    e.inclination_rad=atan2(nlen,W[2]);
    if(!isfinite(ecc) || !isfinite(e.p_m) || !(e.p_m>0)) return TINYELEMENTS_INVALID;
    if(nlen>=1e-12) {
        N[0]=-W[1]/nlen;N[1]=W[0]/nlen;N[2]=0;
        e.raan_rad=wrap(atan2(N[1],N[0]));
        e.arg_periapsis_rad=e.eccentricity==0 ? 0 : angle(N,en,W);
        e.true_anomaly_rad=e.eccentricity==0 ? angle(N,s->position_m,W) : angle(en,s->position_m,W);
    } else {
        double sign=W[2]<0 ? -1.0 : 1.0;
        e.raan_rad=0;
        e.arg_periapsis_rad=e.eccentricity==0 ? 0 : wrap(atan2(sign*en[1],en[0]));
        e.true_anomaly_rad=e.eccentricity==0 ? wrap(atan2(sign*s->position_m[1],s->position_m[0])) : angle(en,s->position_m,W);
    }
    *out=e;return TINYELEMENTS_OK;
}
tinyelements_status tinyelements_from_equinoctial(double mu,const tinyelements_equinoctial *q,tinyelements_state *out)
{
    tinyelements_classical e;double longperi;
    if(!q || !out || !isfinite(q->p_m) || !isfinite(q->f) || !isfinite(q->g) ||
       !isfinite(q->h) || !isfinite(q->k) || !isfinite(q->longitude_rad)) return TINYELEMENTS_INVALID;
    e.p_m=q->p_m;e.eccentricity=hypot(q->f,q->g);e.inclination_rad=2*atan(hypot(q->h,q->k));
    if(TE_PI-e.inclination_rad<1e-8) return TINYELEMENTS_SINGULAR;
    e.raan_rad=(q->h==0 && q->k==0) ? 0 : atan2(q->k,q->h);
    longperi=e.eccentricity<1e-12 ? e.raan_rad : atan2(q->g,q->f);
    e.arg_periapsis_rad=longperi-e.raan_rad;e.true_anomaly_rad=q->longitude_rad-longperi;
    return tinyelements_from_classical(mu,&e,out);
}
tinyelements_status tinyelements_to_equinoctial(double mu,const tinyelements_state *s,tinyelements_equinoctial *out)
{
    tinyelements_classical e;tinyelements_equinoctial q;tinyelements_status status;double t;
    if(!out) return TINYELEMENTS_INVALID;
    status=tinyelements_to_classical(mu,s,&e);if(status!=TINYELEMENTS_OK)return status;
    if(TE_PI-e.inclination_rad<1e-8) return TINYELEMENTS_SINGULAR;
    t=tan(e.inclination_rad/2);q.p_m=e.p_m;
    q.f=e.eccentricity*cos(e.raan_rad+e.arg_periapsis_rad);q.g=e.eccentricity*sin(e.raan_rad+e.arg_periapsis_rad);
    q.h=t*cos(e.raan_rad);q.k=t*sin(e.raan_rad);q.longitude_rad=wrap(e.raan_rad+e.arg_periapsis_rad+e.true_anomaly_rad);
    *out=q;return TINYELEMENTS_OK;
}
tinyelements_status tinyelements_local_frame(double mu,const tinyelements_state *s,tinyelements_frame frame,double out[9])
{
    tinyelements_classical e;tinyelements_status status;double A[3],B[3],W[3],r[9],h[3],n;unsigned j;
    if(!out || (frame!=TINYELEMENTS_RSW && frame!=TINYELEMENTS_NTW && frame!=TINYELEMENTS_PQW)) return TINYELEMENTS_INVALID;
    status=tinyelements_to_classical(mu,s,&e);if(status!=TINYELEMENTS_OK)return status;
    if(frame==TINYELEMENTS_PQW) basis(&e,A,B,W);
    else {
        cross(s->position_m,s->velocity_m_s,h);n=norm(h);for(j=0;j<3;++j)W[j]=h[j]/n;
        if(frame==TINYELEMENTS_RSW) {n=norm(s->position_m);for(j=0;j<3;++j)A[j]=s->position_m[j]/n;cross(W,A,B);}
        else {n=norm(s->velocity_m_s);for(j=0;j<3;++j)B[j]=s->velocity_m_s[j]/n;cross(B,W,A);}
    }
    for(j=0;j<3;++j){r[j]=A[j];r[3+j]=B[j];r[6+j]=W[j];}
    for(j=0;j<9;++j)out[j]=r[j];
    return TINYELEMENTS_OK;
}
