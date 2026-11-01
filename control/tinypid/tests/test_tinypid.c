#include "tinypid/tinypid.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL %d: %s\n",__LINE__,#c);return 1;}}while(0)
#include "data/reference.h"
int main(void){tinypid_state state={{0,0,0}};double command[3],kp[]={.03,.02,.04},kd[]={.01,.015,.02},ki[]={.004,.003,.005},imax[]={.2,.3,.25},tmax[]={.02,.025,.018},maximum=0;unsigned flags;
 for(size_t c=0;c<sizeof reference_cases/sizeof reference_cases[0];++c){CHECK(tinypid_step(&state,reference_cases[c].current,reference_cases[c].desired,reference_cases[c].rate,kp,kd,ki,NULL,.05,imax,tmax,command,&flags)==TINYPID_OK);for(size_t i=0;i<3;++i){double e=fabs(command[i]-reference_cases[c].expected[i]);CHECK(e<3e-14);CHECK(fabs(state.integral[i]-reference_cases[c].integral[i])<3e-14);CHECK(fabs(command[i])<=tmax[i]);CHECK(fabs(state.integral[i])<=imax[i]);maximum=fmax(maximum,e);}}
 {tinypid_state before=state;double saved[]={123,456,789};memcpy(command,saved,sizeof command);CHECK(tinypid_step(&state,reference_cases[0].current,reference_cases[0].desired,reference_cases[0].rate,kp,kd,ki,NULL,-1,imax,tmax,command,&flags)==TINYPID_INVALID_ARGUMENT);CHECK(!memcmp(&state,&before,sizeof state)&&!memcmp(command,saved,sizeof command));}
 {tinypid_state fresh={{0,0,0}};double q[]={cos(.5),0,0,-sin(.5)},desired[]={1,0,0,0},rate[]={0,0,0},large[]={1,1,1},small[]={.00001,.00001,.00001};for(size_t i=0;i<100;++i)CHECK(tinypid_step(&fresh,q,desired,rate,large,kd,ki,NULL,.1,imax,small,command,NULL)==TINYPID_OK);CHECK(fresh.integral[2]==0&&command[2]<0);}
 printf("tinypid: 128 SciPy/NumPy anti-windup cases; max torque error %.3g N m\n",maximum);return 0;}
