/* SPDX-License-Identifier: MIT */
#include "tinyrotation/tinyrotation.h"
#include <float.h>
#include <math.h>
#include <string.h>

static int finite_array(const double *a,size_t n){
    if(!a)return 0;
    for(size_t i=0;i<n;++i)if(!isfinite(a[i]))return 0;
    return 1;
}
tinyrotation_status tinyrotation_normalize(const double q[4],double out[4]){
    double r[4],scale=0,norm=0;
    if(!out||!finite_array(q,4))return TINYROTATION_INVALID_ARGUMENT;
    for(size_t i=0;i<4;++i)scale=fmax(scale,fabs(q[i]));
    if(scale==0)return TINYROTATION_DEGENERATE;
    for(size_t i=0;i<4;++i){r[i]=q[i]/scale;norm+=r[i]*r[i];}
    norm=sqrt(norm);for(size_t i=0;i<4;++i)r[i]/=norm;
    memcpy(out,r,sizeof r);return TINYROTATION_OK;
}
tinyrotation_status tinyrotation_multiply(const double p[4],const double q[4],double out[4]){
    double a[4],b[4],r[4];tinyrotation_status status;
    if(!out)return TINYROTATION_INVALID_ARGUMENT;
    status=tinyrotation_normalize(p,a);if(status!=TINYROTATION_OK)return status;
    status=tinyrotation_normalize(q,b);if(status!=TINYROTATION_OK)return status;
    r[0]=a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3];
    r[1]=a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2];
    r[2]=a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1];
    r[3]=a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0];
    return tinyrotation_normalize(r,out);
}
tinyrotation_status tinyrotation_inverse(const double q[4],double out[4]){
    double r[4];tinyrotation_status status;
    if(!out)return TINYROTATION_INVALID_ARGUMENT;
    status=tinyrotation_normalize(q,r);if(status!=TINYROTATION_OK)return status;
    for(size_t i=1;i<4;++i)r[i]=-r[i];memcpy(out,r,sizeof r);return TINYROTATION_OK;
}
tinyrotation_status tinyrotation_to_dcm(const double q[4],double dcm[9]){
    double a[4],r[9],w,x,y,z;tinyrotation_status status;
    if(!dcm)return TINYROTATION_INVALID_ARGUMENT;
    status=tinyrotation_normalize(q,a);if(status!=TINYROTATION_OK)return status;
    w=a[0];x=a[1];y=a[2];z=a[3];
    r[0]=1-2*(y*y+z*z);r[1]=2*(x*y-w*z);r[2]=2*(x*z+w*y);
    r[3]=2*(x*y+w*z);r[4]=1-2*(x*x+z*z);r[5]=2*(y*z-w*x);
    r[6]=2*(x*z-w*y);r[7]=2*(y*z+w*x);r[8]=1-2*(x*x+y*y);
    memcpy(dcm,r,sizeof r);return TINYROTATION_OK;
}
tinyrotation_status tinyrotation_from_dcm(const double m[9],double q[4]){
    double r[4],s,trace,det;
    if(!q||!finite_array(m,9))return TINYROTATION_INVALID_ARGUMENT;
    for(size_t i=0;i<3;++i)for(size_t j=i;j<3;++j){
        double v=0;for(size_t k=0;k<3;++k)v+=m[k*3+i]*m[k*3+j];
        if(fabs(v-(i==j?1.0:0.0))>1e-10)return TINYROTATION_NOT_A_ROTATION;
    }
    det=m[0]*(m[4]*m[8]-m[5]*m[7])-m[1]*(m[3]*m[8]-m[5]*m[6])+m[2]*(m[3]*m[7]-m[4]*m[6]);
    if(fabs(det-1)>1e-10)return TINYROTATION_NOT_A_ROTATION;
    trace=m[0]+m[4]+m[8];
    if(trace>0){s=2*sqrt(trace+1);r[0]=s/4;r[1]=(m[7]-m[5])/s;r[2]=(m[2]-m[6])/s;r[3]=(m[3]-m[1])/s;}
    else if(m[0]>=m[4]&&m[0]>=m[8]){s=2*sqrt(1+m[0]-m[4]-m[8]);r[0]=(m[7]-m[5])/s;r[1]=s/4;r[2]=(m[1]+m[3])/s;r[3]=(m[2]+m[6])/s;}
    else if(m[4]>=m[8]){s=2*sqrt(1+m[4]-m[0]-m[8]);r[0]=(m[2]-m[6])/s;r[1]=(m[1]+m[3])/s;r[2]=s/4;r[3]=(m[5]+m[7])/s;}
    else{s=2*sqrt(1+m[8]-m[0]-m[4]);r[0]=(m[3]-m[1])/s;r[1]=(m[2]+m[6])/s;r[2]=(m[5]+m[7])/s;r[3]=s/4;}
    if(r[0]<0)for(size_t i=0;i<4;++i)r[i]=-r[i];
    return tinyrotation_normalize(r,q);
}
tinyrotation_status tinyrotation_apply(const double q[4],const double reference[3],double body[3]){
    double m[9],r[3];tinyrotation_status status;
    if(!body||!finite_array(reference,3))return TINYROTATION_INVALID_ARGUMENT;
    status=tinyrotation_to_dcm(q,m);if(status!=TINYROTATION_OK)return status;
    for(size_t i=0;i<3;++i){r[i]=0;for(size_t j=0;j<3;++j)r[i]+=m[i*3+j]*reference[j];}
    if(!finite_array(r,3))return TINYROTATION_ARITHMETIC_ERROR;
    memcpy(body,r,sizeof r);return TINYROTATION_OK;
}
tinyrotation_status tinyrotation_from_axis_angle(const double axis[3],double angle,double q[4]){
    double scale,norm,u[3],r[4],s;
    if(!q||!isfinite(angle)||!finite_array(axis,3))return TINYROTATION_INVALID_ARGUMENT;
    scale=fmax(fabs(axis[0]),fmax(fabs(axis[1]),fabs(axis[2])));if(scale==0)return TINYROTATION_DEGENERATE;
    norm=0;for(size_t i=0;i<3;++i){u[i]=axis[i]/scale;norm+=u[i]*u[i];}norm=sqrt(norm);
    s=sin(angle/2);r[0]=cos(angle/2);for(size_t i=0;i<3;++i)r[i+1]=s*u[i]/norm;
    return tinyrotation_normalize(r,q);
}
tinyrotation_status tinyrotation_to_axis_angle(const double q[4],double axis[3],double *angle){
    double a[4],r[3],theta,norm;tinyrotation_status status;
    if(!axis||!angle)return TINYROTATION_INVALID_ARGUMENT;
    status=tinyrotation_normalize(q,a);if(status!=TINYROTATION_OK)return status;
    if(a[0]<0)for(size_t i=0;i<4;++i)a[i]=-a[i];
    norm=hypot(hypot(a[1],a[2]),a[3]);theta=2*atan2(norm,a[0]);
    if(norm==0){r[0]=1;r[1]=r[2]=0;}else for(size_t i=0;i<3;++i)r[i]=a[i+1]/norm;
    memcpy(axis,r,sizeof r);*angle=theta;return TINYROTATION_OK;
}
tinyrotation_status tinyrotation_from_rotvec(const double vector[3],double q[4]){
    double angle,r[4]={1,0,0,0};
    if(!q||!finite_array(vector,3))return TINYROTATION_INVALID_ARGUMENT;
    angle=hypot(hypot(vector[0],vector[1]),vector[2]);
    if(!isfinite(angle))return TINYROTATION_ARITHMETIC_ERROR;
    if(angle==0){memcpy(q,r,sizeof r);return TINYROTATION_OK;}
    return tinyrotation_from_axis_angle(vector,angle,q);
}
tinyrotation_status tinyrotation_to_rotvec(const double q[4],double vector[3]){
    double axis[3],angle;tinyrotation_status status;
    if(!vector)return TINYROTATION_INVALID_ARGUMENT;
    status=tinyrotation_to_axis_angle(q,axis,&angle);if(status!=TINYROTATION_OK)return status;
    for(size_t i=0;i<3;++i)axis[i]*=angle;memcpy(vector,axis,sizeof axis);return TINYROTATION_OK;
}
tinyrotation_status tinyrotation_from_euler_zyx(double yaw,double pitch,double roll,double q[4]){
    double cy,sy,cp,sp,cr,sr,r[4];
    if(!q||!isfinite(yaw)||!isfinite(pitch)||!isfinite(roll))return TINYROTATION_INVALID_ARGUMENT;
    cy=cos(yaw/2);sy=sin(yaw/2);cp=cos(pitch/2);sp=sin(pitch/2);cr=cos(roll/2);sr=sin(roll/2);
    r[0]=cy*cp*cr+sy*sp*sr;r[1]=cy*cp*sr-sy*sp*cr;r[2]=cy*sp*cr+sy*cp*sr;r[3]=sy*cp*cr-cy*sp*sr;
    return tinyrotation_normalize(r,q);
}
tinyrotation_status tinyrotation_to_euler_zyx(const double q[4],double angles[3]){
    double m[9],r[3],cp;tinyrotation_status status;
    if(!angles)return TINYROTATION_INVALID_ARGUMENT;
    status=tinyrotation_to_dcm(q,m);if(status!=TINYROTATION_OK)return status;
    cp=hypot(m[0],m[3]);r[1]=atan2(-m[6],cp);
    if(cp<=32*DBL_EPSILON){r[0]=atan2(-m[1],m[4]);r[2]=0;}
    else{r[0]=atan2(m[3],m[0]);r[2]=atan2(m[7],m[8]);}
    memcpy(angles,r,sizeof r);return TINYROTATION_OK;
}
tinyrotation_status tinyrotation_slerp(const double from[4],const double to[4],double fraction,double out[4]){
    double a[4],b[4],r[4],dot=0,s0,s1,theta;tinyrotation_status status;
    if(!out||!isfinite(fraction)||fraction<0||fraction>1)return TINYROTATION_INVALID_ARGUMENT;
    status=tinyrotation_normalize(from,a);if(status!=TINYROTATION_OK)return status;
    status=tinyrotation_normalize(to,b);if(status!=TINYROTATION_OK)return status;
    for(size_t i=0;i<4;++i)dot+=a[i]*b[i];
    if(dot<0){dot=-dot;for(size_t i=0;i<4;++i)b[i]=-b[i];}
    dot=fmin(dot,1);
    if(dot>1-1e-12){s0=1-fraction;s1=fraction;}
    else{theta=acos(dot);s0=sin((1-fraction)*theta)/sin(theta);s1=sin(fraction*theta)/sin(theta);}
    for(size_t i=0;i<4;++i)r[i]=s0*a[i]+s1*b[i];return tinyrotation_normalize(r,out);
}
