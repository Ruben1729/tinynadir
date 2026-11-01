#include "tinydoubler/tinydoubler.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL %d: %s\n",__LINE__,#c);return 1;}}while(0)
#include "data/reference.h"
int main(void){tinydoubler_result out;double max_r=0,max_v=0;unsigned validated_vallado=0;
 for(size_t c=0;c<sizeof reference_cases/sizeof reference_cases[0];++c){CHECK(tinydoubler_estimate(reference_cases[c].mu,reference_cases[c].t1,reference_cases[c].t3,reference_cases[c].los,reference_cases[c].sites,reference_cases[c].radius1,reference_cases[c].radius2,1e-9,64,&out)==TINYDOUBLER_OK);CHECK(out.iterations<=64);CHECK(fabs(out.time_residual_s[0])<=1e-9&&fabs(out.time_residual_s[1])<=1e-9);double er=0,ev=0;for(size_t i=0;i<3;++i){double dr=out.position_m[i]-reference_cases[c].expected_position[i],dv=out.velocity_m_s[i]-reference_cases[c].expected_velocity[i];er+=dr*dr;ev+=dv*dv;CHECK(out.ranges_m[i]>0);if(reference_cases[c].vallado_valid){CHECK(fabs(out.position_m[i]-reference_cases[c].vallado_position[i])<1e-3);CHECK(fabs(out.velocity_m_s[i]-reference_cases[c].vallado_velocity[i])<1e-4);}}er=sqrt(er);ev=sqrt(ev);CHECK(er<1e-3&&ev<1e-4);max_r=fmax(max_r,er);max_v=fmax(max_v,ev);validated_vallado+=reference_cases[c].vallado_valid;}
 CHECK(validated_vallado==37);
 {tinydoubler_result before=out;CHECK(tinydoubler_estimate(reference_cases[0].mu,reference_cases[0].t1,reference_cases[0].t3,reference_cases[0].los,reference_cases[0].sites,reference_cases[0].radius1,reference_cases[0].radius2,1e-10,1,&out)==TINYDOUBLER_NO_CONVERGENCE);CHECK(!memcmp(&out,&before,sizeof out));
 double bad_los[9];memcpy(bad_los,reference_cases[0].los,sizeof bad_los);bad_los[4]=NAN;CHECK(tinydoubler_estimate(reference_cases[0].mu,reference_cases[0].t1,reference_cases[0].t3,bad_los,reference_cases[0].sites,reference_cases[0].radius1,reference_cases[0].radius2,1e-9,64,&out)==TINYDOUBLER_INVALID_ARGUMENT);CHECK(!memcmp(&out,&before,sizeof out));}
 {tinydoubler_result before=out;double parallel[]={1,0,0,1,0,0,1,0,0};CHECK(tinydoubler_estimate(reference_cases[0].mu,-20,20,parallel,reference_cases[0].sites,1,1,1e-9,64,&out)==TINYDOUBLER_GEOMETRY);CHECK(!memcmp(&out,&before,sizeof out));CHECK(tinydoubler_estimate(reference_cases[0].mu,0,20,reference_cases[0].los,reference_cases[0].sites,7e6,7e6,1e-9,64,&out)==TINYDOUBLER_INVALID_ARGUMENT);CHECK(tinydoubler_estimate(reference_cases[0].mu,-20,20,reference_cases[0].los,reference_cases[0].sites,7e6,7e6,1e-9,65,&out)==TINYDOUBLER_INVALID_ARGUMENT);CHECK(tinydoubler_estimate(reference_cases[0].mu,-20,20,NULL,reference_cases[0].sites,7e6,7e6,1e-9,64,&out)==TINYDOUBLER_INVALID_ARGUMENT);}
 printf("tinydoubler: 48 SciPy conic fits / 37 validated Vallado cases; max %.3g m / %.3g m/s\n",max_r,max_v);return 0;}
