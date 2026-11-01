/* SPDX-License-Identifier: MIT */
#include "tinydoubler/tinydoubler.h"
#include "tinylinalg/tinylinalg.h"
#include "tinygibbs/tinygibbs.h"
#include <math.h>
#include <stddef.h>
#include <string.h>
#define TD_PI 3.14159265358979323846264338327950288

static int finite_array(const double *a,size_t n){if(!a)return 0;for(size_t i=0;i<n;++i)if(!isfinite(a[i]))return 0;return 1;}
static double advance(double angle){angle=fmod(angle,2*TD_PI);if(angle<=0)angle+=2*TD_PI;return angle;}
static int evaluate(double mu,double t1,double t3,const double los[9],const double sites[9],
    double radius1,double radius2,tinydoubler_result *result){
    tinydoubler_result next={0};double position[9],normal[3],cross[3],radii[2]={radius1,radius2};
    double h[3],hhat[3],evec[3],r2,alpha,a,e,mean_motion,mean_anomaly[3],dm[2];tinygibbs_result velocity;
    for(size_t i=0;i<2;++i){
        double projection=tinylinalg_dot3(los+i*3,sites+i*3);
        double discriminant=projection*projection+radii[i]*radii[i]-tinylinalg_dot3(sites+i*3,sites+i*3);
        if(!isfinite(discriminant)||discriminant<=0)return 0;
        next.ranges_m[i]=-projection+sqrt(discriminant);
        if(next.ranges_m[i]<=0||!isfinite(next.ranges_m[i]))return 0;
        for(size_t j=0;j<3;++j)position[i*3+j]=sites[i*3+j]+next.ranges_m[i]*los[i*3+j];
    }
    if(tinylinalg_cross3(position,position+3,cross)!=TINYLINALG_OK||tinylinalg_unit3(cross,normal)!=TINYLINALG_OK)return 0;
    double denominator=tinylinalg_dot3(normal,los+6);
    if(fabs(denominator)<=1e-12)return 0;
    next.ranges_m[2]=-tinylinalg_dot3(normal,sites+6)/denominator;
    if(!isfinite(next.ranges_m[2])||next.ranges_m[2]<=0)return 0;
    for(size_t j=0;j<3;++j)position[6+j]=sites[6+j]+next.ranges_m[2]*los[6+j];
    if(!finite_array(position,9)||tinygibbs_estimate(mu,position,position+3,position+6,1e-6,&velocity)!=TINYGIBBS_OK)return 0;
    r2=tinylinalg_norm3(position+3);alpha=2/r2-tinylinalg_dot3(velocity.velocity_m_s,velocity.velocity_m_s)/mu;
    if(!isfinite(alpha)||alpha<=0)return 0;
    a=1/alpha;mean_motion=sqrt(mu/a)/a;if(!isfinite(mean_motion)||mean_motion<=0)return 0;
    if(tinylinalg_cross3(position+3,velocity.velocity_m_s,h)!=TINYLINALG_OK||tinylinalg_unit3(h,hhat)!=TINYLINALG_OK||tinylinalg_cross3(velocity.velocity_m_s,h,cross)!=TINYLINALG_OK)return 0;
    for(size_t i=0;i<3;++i)evec[i]=cross[i]/mu-position[3+i]/r2;
    e=tinylinalg_norm3(evec);if(!isfinite(e)||e>=.95)return 0;
    for(size_t i=0;i<2;++i){
        if(tinylinalg_cross3(position+i*3,position+(i+1)*3,cross)!=TINYLINALG_OK||tinylinalg_dot3(cross,hhat)<=0)return 0;
    }
    if(e>1e-12){
        double p[3],q[3];for(size_t i=0;i<3;++i)p[i]=evec[i]/e;
        if(tinylinalg_cross3(hhat,p,q)!=TINYLINALG_OK||tinylinalg_unit3(q,q)!=TINYLINALG_OK||tinylinalg_cross3(q,hhat,p)!=TINYLINALG_OK)return 0;
        for(size_t i=0;i<3;++i){double cosine=tinylinalg_dot3(position+i*3,p)/a+e,sine=tinylinalg_dot3(position+i*3,q)/(a*sqrt(1-e*e));double eccentric=atan2(sine,cosine);mean_anomaly[i]=eccentric-e*sin(eccentric);}
        dm[0]=advance(mean_anomaly[1]-mean_anomaly[0]);dm[1]=advance(mean_anomaly[2]-mean_anomaly[1]);
    }else{
        for(size_t i=0;i<2;++i){if(tinylinalg_cross3(position+i*3,position+(i+1)*3,cross)!=TINYLINALG_OK)return 0;dm[i]=advance(atan2(tinylinalg_dot3(hhat,cross),tinylinalg_dot3(position+i*3,position+(i+1)*3)));}
    }
    if(dm[0]>TD_PI||dm[1]>TD_PI)return 0;
    next.time_residual_s[0]=dm[0]/mean_motion+t1;next.time_residual_s[1]=dm[1]/mean_motion-t3;
    if(!finite_array(next.time_residual_s,2))return 0;
    memcpy(next.position_m,position+3,sizeof next.position_m);memcpy(next.velocity_m_s,velocity.velocity_m_s,sizeof next.velocity_m_s);next.radii_m[0]=radius1;next.radii_m[1]=radius2;*result=next;return 1;
}
tinydoubler_status tinydoubler_estimate(double mu,double t1,double t3,const double los[9],const double sites[9],
    double radius1,double radius2,double tolerance,unsigned maximum_iterations,tinydoubler_result *result){
    double directions[9],scale=0,x[2];tinydoubler_result current;
    if(!result||!isfinite(mu)||mu<=0||!isfinite(t1)||t1>=0||fabs(t1)>1e6||!isfinite(t3)||t3<=0||t3>1e6||!isfinite(radius1)||radius1<=0||!isfinite(radius2)||radius2<=0||!isfinite(tolerance)||tolerance<1e-10||tolerance>1||maximum_iterations<1||maximum_iterations>64||!finite_array(los,9)||!finite_array(sites,9))return TINYDOUBLER_INVALID_ARGUMENT;
    for(size_t i=0;i<3;++i){if(tinylinalg_unit3(los+3*i,directions+3*i)!=TINYLINALG_OK)return TINYDOUBLER_GEOMETRY;scale=fmax(scale,tinylinalg_norm3(sites+3*i));}
    if(!isfinite(scale)||scale==0)return TINYDOUBLER_GEOMETRY;
    x[0]=radius1/scale;x[1]=radius2/scale;
    if(!evaluate(mu,t1,t3,directions,sites,radius1,radius2,&current))return TINYDOUBLER_GEOMETRY;
    for(unsigned iteration=0;iteration<=maximum_iterations;++iteration){
        double residual_norm=hypot(current.time_residual_s[0],current.time_residual_s[1]),jacobian[4],step[2],determinant;
        if(fmax(fabs(current.time_residual_s[0]),fabs(current.time_residual_s[1]))<=tolerance){current.iterations=iteration;*result=current;return TINYDOUBLER_OK;}
        if(iteration==maximum_iterations)break;
        for(size_t j=0;j<2;++j){
            double h=1e-5*fmax(1,fabs(x[j])),plus[2]={x[0],x[1]},minus[2]={x[0],x[1]};tinydoubler_result p,m;
            plus[j]+=h;minus[j]-=h;
            int pv=evaluate(mu,t1,t3,directions,sites,plus[0]*scale,plus[1]*scale,&p);
            int mv=minus[j]>0&&evaluate(mu,t1,t3,directions,sites,minus[0]*scale,minus[1]*scale,&m);
            if(!pv&&!mv)return TINYDOUBLER_NO_CONVERGENCE;
            for(size_t i=0;i<2;++i)jacobian[i*2+j]=pv&&mv?(p.time_residual_s[i]-m.time_residual_s[i])/(2*h):pv?(p.time_residual_s[i]-current.time_residual_s[i])/h:(current.time_residual_s[i]-m.time_residual_s[i])/h;
        }
        determinant=jacobian[0]*jacobian[3]-jacobian[1]*jacobian[2];
        if(!isfinite(determinant)||fabs(determinant)<=1e-12*hypot(jacobian[0],jacobian[1])*hypot(jacobian[2],jacobian[3]))return TINYDOUBLER_GEOMETRY;
        step[0]=(jacobian[3]*current.time_residual_s[0]-jacobian[1]*current.time_residual_s[1])/determinant;
        step[1]=(-jacobian[2]*current.time_residual_s[0]+jacobian[0]*current.time_residual_s[1])/determinant;
        if(!finite_array(step,2))return TINYDOUBLER_NO_CONVERGENCE;
        double fraction=1,relative=hypot(step[0],step[1])/hypot(x[0],x[1]);if(relative>.25)fraction=.25/relative;
        int improved=0;
        for(size_t search=0;search<12;++search){
            double trial[2]={x[0]-fraction*step[0],x[1]-fraction*step[1]};tinydoubler_result candidate;
            if(trial[0]>0&&trial[1]>0&&evaluate(mu,t1,t3,directions,sites,trial[0]*scale,trial[1]*scale,&candidate)&&hypot(candidate.time_residual_s[0],candidate.time_residual_s[1])<residual_norm){x[0]=trial[0];x[1]=trial[1];current=candidate;improved=1;break;}
            fraction*=.5;
        }
        if(!improved)return TINYDOUBLER_NO_CONVERGENCE;
    }
    return TINYDOUBLER_NO_CONVERGENCE;
}
