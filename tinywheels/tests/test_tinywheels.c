#include "tinywheels/tinywheels.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL %d: %s\n",__LINE__,#c);return 1;}}while(0)
#include "data/reference.h"
int main(void){tinywheels_workspace work;double command[8],residual[3],maximum=0;
 for(size_t c=0;c<sizeof reference_cases/sizeof reference_cases[0];++c){CHECK(tinywheels_allocate(reference_cases[c].count,reference_cases[c].axes,reference_cases[c].torque,reference_cases[c].speed,reference_cases[c].tlimit,reference_cases[c].slimit,reference_cases[c].enabled,command,residual,&work)==TINYWHEELS_OK);for(size_t i=0;i<reference_cases[c].count;++i){double e=fabs(command[i]-reference_cases[c].expected[i]);CHECK(e<3e-12);CHECK(fabs(command[i])<=reference_cases[c].tlimit[i]);maximum=fmax(maximum,e);}for(size_t i=0;i<3;++i)CHECK(fabs(residual[i]-reference_cases[c].residual[i])<3e-12);}
 {double axes[]={1,0,0,0,1,0,0,0,1},torque[]={.1,0,0},speed[]={0,0,0},limit[]={1,1,1},slimit[]={10,10,10};CHECK(tinywheels_allocate(3,axes,torque,speed,limit,slimit,NULL,command,residual,&work)==TINYWHEELS_OK);CHECK(fabs(command[0]+.1)<1e-14&&fabs(residual[0])<1e-14);speed[0]=-10;CHECK(tinywheels_allocate(3,axes,torque,speed,limit,slimit,NULL,command,residual,&work)==TINYWHEELS_OK);CHECK(command[0]==0&&fabs(residual[0]-.1)<1e-14);}
 {double a[]={1,0,0},t[]={1,0,0},s[]={0},l[]={0},sl[]={10},saved[]={123,456,789};memcpy(residual,saved,sizeof residual);command[0]=42;CHECK(tinywheels_allocate(1,a,t,s,l,sl,NULL,command,residual,&work)==TINYWHEELS_NO_AUTHORITY);CHECK(command[0]==42&&!memcmp(residual,saved,sizeof residual));CHECK(tinywheels_allocate(9,a,t,s,l,sl,NULL,command,residual,&work)==TINYWHEELS_INVALID_ARGUMENT);}
 printf("tinywheels: 96 NumPy SVD allocation cases; max wheel torque error %.3g N m\n",maximum);return 0;}
