/* SPDX-License-Identifier: MIT */
#include "tinyrotation/tinyrotation.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL %s:%d: %s\n",__FILE__,__LINE__,#c);return 1;}}while(0)
#include "data/reference.h"
static double qdistance(const double a[4],const double b[4]){
    double minus=0,plus=0;for(size_t i=0;i<4;++i){minus+=(a[i]-b[i])*(a[i]-b[i]);plus+=(a[i]+b[i])*(a[i]+b[i]);}return sqrt(fmin(minus,plus));
}
static int close_array(const double *a,const double *b,size_t n,double tol){for(size_t i=0;i<n;++i)if(fabs(a[i]-b[i])>tol*fmax(1,fabs(b[i])))return 0;return 1;}
int main(void){
    double q[4],m[9],r[3],axis[3],angle,angles[3],max_matrix=0,max_quaternion=0;
    for(size_t c=0;c<sizeof reference_cases/sizeof reference_cases[0];++c){
        const double *expected=reference_cases[c].q;
        CHECK(tinyrotation_to_dcm(expected,m)==TINYROTATION_OK);
        CHECK(close_array(m,reference_cases[c].matrix,9,3e-14));
        for(size_t i=0;i<9;++i)max_matrix=fmax(max_matrix,fabs(m[i]-reference_cases[c].matrix[i]));
        CHECK(tinyrotation_from_dcm(m,q)==TINYROTATION_OK);CHECK(qdistance(q,expected)<1e-14);
        CHECK(tinyrotation_apply(expected,reference_cases[c].vector,r)==TINYROTATION_OK);CHECK(close_array(r,reference_cases[c].rotated,3,3e-14));
        CHECK(tinyrotation_to_rotvec(expected,r)==TINYROTATION_OK);CHECK(close_array(r,reference_cases[c].rotvec,3,3e-14));
        CHECK(tinyrotation_from_rotvec(r,q)==TINYROTATION_OK);CHECK(qdistance(q,expected)<1e-14);
        CHECK(tinyrotation_to_axis_angle(expected,axis,&angle)==TINYROTATION_OK);
        CHECK(angle>=0&&angle<=3.141592653589794);
        CHECK(tinyrotation_from_axis_angle(axis,angle,q)==TINYROTATION_OK);CHECK(qdistance(q,expected)<1e-14);
        CHECK(tinyrotation_multiply(expected,reference_cases[c].q2,q)==TINYROTATION_OK);
        CHECK(qdistance(q,reference_cases[c].composition)<1e-14);
        CHECK(tinyrotation_slerp(expected,reference_cases[c].q2,reference_cases[c].fraction,q)==TINYROTATION_OK);
        CHECK(qdistance(q,reference_cases[c].slerp)<1e-14);max_quaternion=fmax(max_quaternion,qdistance(q,reference_cases[c].slerp));
        CHECK(tinyrotation_to_euler_zyx(expected,angles)==TINYROTATION_OK);
        CHECK(tinyrotation_from_euler_zyx(angles[0],angles[1],angles[2],q)==TINYROTATION_OK);CHECK(qdistance(q,expected)<3e-13);
        CHECK(tinyrotation_inverse(expected,q)==TINYROTATION_OK);
        CHECK(tinyrotation_multiply(expected,q,q)==TINYROTATION_OK);CHECK(fabs(fabs(q[0])-1)<1e-14);
        CHECK(fabs(q[1])+fabs(q[2])+fabs(q[3])<1e-14);
    }
    {double zero[4]={0},bad[4]={1,NAN,0,0},out[4]={123,456,789,1},saved[4];memcpy(saved,out,sizeof out);
      CHECK(tinyrotation_normalize(zero,out)==TINYROTATION_DEGENERATE);CHECK(tinyrotation_normalize(bad,out)==TINYROTATION_INVALID_ARGUMENT);
      CHECK(tinyrotation_slerp(bad,bad,0.5,out)==TINYROTATION_INVALID_ARGUMENT);CHECK(!memcmp(saved,out,sizeof out));
      CHECK(tinyrotation_slerp(saved,saved,-0.1,out)==TINYROTATION_INVALID_ARGUMENT);CHECK(tinyrotation_normalize(NULL,out)==TINYROTATION_INVALID_ARGUMENT);}
    {double large[4]={DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX};CHECK(tinyrotation_normalize(large,q)==TINYROTATION_OK);for(size_t i=0;i<4;++i)CHECK(q[i]==0.5);}
    {double reflection[9]={1,0,0,0,1,0,0,0,-1},bad[9]={1,0.1,0,0,1,0,0,0,1};CHECK(tinyrotation_from_dcm(reflection,q)==TINYROTATION_NOT_A_ROTATION);CHECK(tinyrotation_from_dcm(bad,q)==TINYROTATION_NOT_A_ROTATION);}
    {double q1[4]={0,1,0,0},q2[4]={0,-1,0,0};CHECK(tinyrotation_slerp(q1,q2,.25,q)==TINYROTATION_OK);CHECK(qdistance(q,q1)<1e-15);}
    {const double pi=3.14159265358979323846;for(int sign=-1;sign<=1;sign+=2){double original[4];CHECK(tinyrotation_from_euler_zyx(.7,sign*pi/2,-.9,original)==TINYROTATION_OK);CHECK(tinyrotation_to_euler_zyx(original,angles)==TINYROTATION_OK);CHECK(angles[2]==0);CHECK(tinyrotation_from_euler_zyx(angles[0],angles[1],angles[2],q)==TINYROTATION_OK);CHECK(qdistance(q,original)<1e-14);}}
    printf("tinyrotation: 96 SciPy cases; max DCM error %.3g, SLERP quaternion error %.3g\n",max_matrix,max_quaternion);return 0;
}
