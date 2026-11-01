#include "tinypd/tinypd.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL %d: %s\n",__LINE__,#c);return 1;}}while(0)
#include "data/reference.h"
int main(void){double out[3],maximum=0;
 for(size_t c=0;c<sizeof reference_cases/sizeof reference_cases[0];++c){CHECK(tinypd_torque(reference_cases[c].current,reference_cases[c].desired,reference_cases[c].rate,reference_cases[c].kp,reference_cases[c].kd,reference_cases[c].ff,out)==TINYPD_OK);for(size_t i=0;i<3;++i){double e=fabs(out[i]-reference_cases[c].expected[i]);CHECK(e<3e-14);maximum=fmax(maximum,e);}}
 {double current[]={cos(.05),0,0,-sin(.05)},desired[]={1,0,0,0},rate[]={0,0,0},kp[]={1,1,1},kd[]={0,0,0};CHECK(tinypd_torque(current,desired,rate,kp,kd,NULL,out)==TINYPD_OK);CHECK(out[2]<0);current[0]=-current[0];current[3]=-current[3];double other[3];CHECK(tinypd_torque(current,desired,rate,kp,kd,NULL,other)==TINYPD_OK);CHECK(fabs(other[2]-out[2])<1e-15);}
 {double q[]={1,0,0,0},v[]={0,0,0},g[]={-1,1,1},saved[]={123,456,789};memcpy(out,saved,sizeof out);CHECK(tinypd_torque(q,q,v,g,v,NULL,out)==TINYPD_INVALID_ARGUMENT);CHECK(!memcmp(out,saved,sizeof out));CHECK(tinypd_torque(NULL,q,v,v,v,NULL,out)==TINYPD_INVALID_ARGUMENT);}
 printf("tinypd: 64 SciPy quaternion feedback cases; max torque error %.3g N m\n",maximum);return 0;}
