#include "tinylaplace/tinylaplace.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL %d: %s\n",__LINE__,#c);return 1;}}while(0)
#include "data/reference.h"
int main(void){tinylaplace_workspace work;tinylaplace_result out;double max_r=0,max_v=0;
 for(size_t c=0;c<sizeof reference_cases/sizeof reference_cases[0];++c){CHECK(tinylaplace_three(reference_cases[c].mu,reference_cases[c].spacing,reference_cases[c].los,reference_cases[c].site,reference_cases[c].observer_velocity,reference_cases[c].observer_acceleration,0,&out,&work)==TINYLAPLACE_OK);CHECK(out.admissible_roots>0&&out.range_m>0);double er=0,ev=0;for(size_t i=0;i<3;++i){double dr=out.position_m[i]-reference_cases[c].expected_position[i],dv=out.velocity_m_s[i]-reference_cases[c].expected_velocity[i];er+=dr*dr;ev+=dv*dv;}er=sqrt(er);ev=sqrt(ev);CHECK(er<1e-4&&ev<1e-5);max_r=fmax(max_r,er);max_v=fmax(max_v,ev);}
 {tinylaplace_result before=out;double l[]={1,0,0},zero[]={0,0,0};CHECK(tinylaplace_estimate(reference_cases[0].mu,l,zero,zero,reference_cases[0].site,reference_cases[0].observer_velocity,reference_cases[0].observer_acceleration,0,&out,&work)==TINYLAPLACE_GEOMETRY);CHECK(!memcmp(&out,&before,sizeof out));CHECK(tinylaplace_three(reference_cases[0].mu,0,reference_cases[0].los,reference_cases[0].site,zero,zero,0,&out,&work)==TINYLAPLACE_INVALID_ARGUMENT);CHECK(tinylaplace_three(reference_cases[0].mu,20,NULL,reference_cases[0].site,zero,zero,0,&out,&work)==TINYLAPLACE_INVALID_ARGUMENT);}
 printf("tinylaplace: 24 Vallado curvature oracle cases; max %.3g m / %.3g m/s\n",max_r,max_v);return 0;}
