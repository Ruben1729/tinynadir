#include "tinyrigidbody/tinyrigidbody.h"
#include "tinyrotation/tinyrotation.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL %d: %s\n",__LINE__,#c);return 1;}}while(0)
#include "data/reference.h"
#include "data/trajectory.h"
static int rhs(const double y[7],double dy[7]){const double j[]={.02,0,0,0,.03,0,0,0,.05},zero[]={0,0,0};tinyrigidbody_workspace work;return tinyrigidbody_derivative(y,y+4,j,zero,NULL,dy,dy+4,&work)!=TINYRIGIDBODY_OK;}
static int step(double y[7],double dt){double k[4][7],tmp[7];if(rhs(y,k[0]))return 1;for(size_t i=0;i<7;++i)tmp[i]=y[i]+dt*k[0][i]/2;if(rhs(tmp,k[1]))return 1;for(size_t i=0;i<7;++i)tmp[i]=y[i]+dt*k[1][i]/2;if(rhs(tmp,k[2]))return 1;for(size_t i=0;i<7;++i)tmp[i]=y[i]+dt*k[2][i];if(rhs(tmp,k[3]))return 1;for(size_t i=0;i<7;++i)y[i]+=dt*(k[0][i]+2*k[1][i]+2*k[2][i]+k[3][i])/6;return tinyrotation_normalize(y,y)!=TINYROTATION_OK;}
int main(void){tinyrigidbody_workspace work;double dq[4],dw[3],maximum=0;
 for(size_t c=0;c<sizeof reference_cases/sizeof reference_cases[0];++c){CHECK(tinyrigidbody_derivative(reference_cases[c].q,reference_cases[c].omega,reference_cases[c].inertia,reference_cases[c].torque,reference_cases[c].momentum,dq,dw,&work)==TINYRIGIDBODY_OK);for(size_t i=0;i<4;++i)CHECK(fabs(dq[i]-reference_cases[c].qdot[i])<3e-10);for(size_t i=0;i<3;++i){double e=fabs(dw[i]-reference_cases[c].omegadot[i]);CHECK(e<3e-13);maximum=fmax(maximum,e);}}
 {double y[]={1,0,0,0,.3,.2,.1};size_t steps=0;for(size_t c=0;c<sizeof trajectory_cases/sizeof trajectory_cases[0];++c){size_t end=(size_t)(trajectory_cases[c].t*1000+.5);while(steps<end){CHECK(step(y,.001)==0);++steps;}for(size_t i=0;i<7;++i)CHECK(fabs(y[i]-trajectory_cases[c].state[i])<3e-11);double energy=.5*(.02*y[4]*y[4]+.03*y[5]*y[5]+.05*y[6]*y[6]);CHECK(fabs(energy-.00175)<3e-14);double h[]={.02*y[4],.03*y[5],.05*y[6]},inverse[4],inertial[3];CHECK(tinyrotation_inverse(y,inverse)==TINYROTATION_OK);CHECK(tinyrotation_apply(inverse,h,inertial)==TINYROTATION_OK);CHECK(fabs(inertial[0]-.006)<3e-13&&fabs(inertial[1]-.006)<3e-13&&fabs(inertial[2]-.005)<3e-13);}}
 {double q[]={1,0,0,0},w[]={0,0,0},j[]={1,0,0,0,-1,0,0,0,1},t[]={0,0,0},saved[]={123,456,789,42};memcpy(dq,saved,sizeof dq);dw[0]=99;CHECK(tinyrigidbody_derivative(q,w,j,t,NULL,dq,dw,&work)==TINYRIGIDBODY_INVALID_INERTIA);CHECK(!memcmp(dq,saved,sizeof dq)&&dw[0]==99);}
 printf("tinyrigidbody: 64 NumPy Euler derivatives + SciPy DOP853 torque-free trajectory; max angular-acceleration error %.3g\n",maximum);return 0;}
