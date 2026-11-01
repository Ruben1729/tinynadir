/* SPDX-License-Identifier: MIT */
#include "tinygauss/tinygauss.h"
#include "tinylinalg/tinylinalg.h"
#include "tinygibbs/tinygibbs.h"
#include <math.h>
#include <stddef.h>
#include <string.h>
static int finite_array(const double *v,size_t n){if(!v)return 0;for(size_t i=0;i<n;++i)if(!isfinite(v[i]))return 0;return 1;}
static double determinant(const double m[9]){return m[0]*(m[4]*m[8]-m[5]*m[7])-m[1]*(m[3]*m[8]-m[5]*m[6])+m[2]*(m[3]*m[7]-m[4]*m[6]);}
tinygauss_status tinygauss_estimate(double mu,double t1,double t3,const double los[9],const double observers[9],
    double hint,tinygauss_result *result,tinygauss_workspace *work){
    double directions[9],sites[9],matrix[9],n[9],scale=0,a1,a3,b1,b3,A,B,projection,mu_scaled,coefficients[9]={1,0,0,0,0,0,0,0,0};
    tinygauss_result best={0};double best_score=HUGE_VAL;unsigned candidates=0;
    if(!result||!work||!isfinite(mu)||mu<=0||!isfinite(t1)||t1>=0||fabs(t1)>1e6||!isfinite(t3)||t3<=0||t3>1e6||!isfinite(hint)||hint<0||!finite_array(los,9)||!finite_array(observers,9))return TINYGAUSS_INVALID_ARGUMENT;
    for(size_t i=0;i<3;++i){if(tinylinalg_unit3(los+3*i,directions+3*i)!=TINYLINALG_OK)return TINYGAUSS_GEOMETRY;scale=fmax(scale,tinylinalg_norm3(observers+3*i));}
    if(!isfinite(scale)||scale==0)return TINYGAUSS_GEOMETRY;
    for(size_t i=0;i<9;++i)sites[i]=observers[i]/scale;
    for(size_t i=0;i<3;++i)for(size_t j=0;j<3;++j)matrix[i*3+j]=directions[j*3+i];
    if(fabs(determinant(matrix))<=1e-12)return TINYGAUSS_GEOMETRY;
    for(size_t j=0;j<3;++j){double rhs[3],solution[3];for(size_t i=0;i<3;++i)rhs[i]=sites[j*3+i];if(tinylinalg_solve(3,matrix,rhs,solution,work->linear,12)!=TINYLINALG_OK)return TINYGAUSS_GEOMETRY;for(size_t i=0;i<3;++i)n[i*3+j]=solution[i];}
    a1=t3/(t3-t1);a3=-t1/(t3-t1);b1=a1*((t3-t1)*(t3-t1)-t3*t3)/6;b3=a3*((t3-t1)*(t3-t1)-t1*t1)/6;
    mu_scaled=((mu/scale)/scale)/scale;A=a1*n[3]+a3*n[5]-n[4];B=mu_scaled*(b1*n[3]+b3*n[5]);projection=tinylinalg_dot3(directions+3,sites+3);
    coefficients[2]=-(tinylinalg_dot3(sites+3,sites+3)+A*A+2*A*projection);coefficients[5]=-2*B*(A+projection);coefficients[8]=-B*B;
    if(!finite_array(coefficients,9))return TINYGAUSS_ARITHMETIC_ERROR;
    if(tinypolynomial_roots(8,coefficients,work->roots,work->scratch_roots,16)!=TINYPOLYNOMIAL_OK)return TINYGAUSS_NO_CONVERGENCE;
    for(size_t root=0;root<8;++root){
        double radius=work->roots[root].real,ranges[3],rhs[3],positions[9],c1,c3,u,score;tinygibbs_result velocity;tinygauss_result candidate={0};
        if(radius<=0||fabs(work->roots[root].imaginary)>1e-8*fmax(1,fabs(radius)))continue;
        u=mu_scaled/radius/radius/radius;c1=a1+b1*u;c3=a3+b3*u;
        for(size_t i=0;i<3;++i){matrix[i*3]=c1*directions[i];matrix[i*3+1]=-directions[3+i];matrix[i*3+2]=c3*directions[6+i];rhs[i]=sites[3+i]-c1*sites[i]-c3*sites[6+i];}
        if(tinylinalg_solve(3,matrix,rhs,ranges,work->linear,12)!=TINYLINALG_OK||ranges[0]<=0||ranges[1]<=0||ranges[2]<=0)continue;
        for(size_t i=0;i<3;++i){candidate.ranges_m[i]=ranges[i]*scale;for(size_t j=0;j<3;++j)positions[i*3+j]=(sites[i*3+j]+ranges[i]*directions[i*3+j])*scale;}
        if(!finite_array(positions,9)||tinygibbs_estimate(mu,positions,positions+3,positions+6,1e-5,&velocity)!=TINYGIBBS_OK)continue;
        memcpy(candidate.position_m,positions+3,sizeof candidate.position_m);memcpy(candidate.velocity_m_s,velocity.velocity_m_s,sizeof candidate.velocity_m_s);candidate.selected_radius_m=radius*scale;
        if(!isfinite(candidate.selected_radius_m))continue;
        ++candidates;score=hint>0?fabs(radius-hint/scale):-radius;
        if(score<best_score){best=candidate;best_score=score;}
    }
    if(!candidates)return TINYGAUSS_NO_SOLUTION;
    best.admissible_roots=candidates;*result=best;return TINYGAUSS_OK;
}
