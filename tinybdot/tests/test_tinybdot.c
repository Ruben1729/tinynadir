#include "tinybdot/tinybdot.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL %d: %s\n",__LINE__,#c);return 1;}}while(0)
#include "data/reference.h"
int main(void){tinybdot_state state={0};double command[3],limit[]={.15,.1,.2},maximum=0;unsigned flags;
 for(size_t c=0;c<sizeof reference_cases/sizeof reference_cases[0];++c){CHECK(tinybdot_step(&state,reference_cases[c].field,.02,12000,.12,limit,command,&flags)==TINYBDOT_OK);for(size_t i=0;i<3;++i){double e=fabs(command[i]-reference_cases[c].expected[i]);CHECK(e<3e-13);CHECK(fabs(state.filtered_derivative[i]-reference_cases[c].filtered[i])<3e-18);CHECK(fabs(command[i])<=limit[i]);maximum=fmax(maximum,e);}}
 {tinybdot_state before=state;double saved[]={123,456,789};memcpy(command,saved,sizeof command);CHECK(tinybdot_step(&state,reference_cases[0].field,0,12000,.12,limit,command,&flags)==TINYBDOT_INVALID_ARGUMENT);CHECK(!memcmp(&state,&before,sizeof state)&&!memcmp(command,saved,sizeof command));}
 {tinybdot_state fresh={0};double b1[]={1e-5,2e-5,3e-5},b2[]={2e-5,2e-5,3e-5};CHECK(tinybdot_step(&fresh,b1,1,1,0,limit,command,NULL)==TINYBDOT_OK);CHECK(command[0]==0);CHECK(tinybdot_step(&fresh,b2,1,1,0,limit,command,NULL)==TINYBDOT_OK);CHECK(fabs(command[0]+1e-5)<1e-20);}
 printf("tinybdot: 128 SciPy lfilter samples; max dipole error %.3g A m^2\n",maximum);return 0;}
