#include "tinygauss/tinygauss.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL %d: %s\n",__LINE__,#c);return 1;}}while(0)
#include "data/reference.h"
int main(void){tinygauss_workspace work;tinygauss_result out;double max_r=0,max_v=0;
 for(size_t c=0;c<sizeof reference_cases/sizeof reference_cases[0];++c){CHECK(tinygauss_estimate(reference_cases[c].mu,reference_cases[c].t1,reference_cases[c].t3,reference_cases[c].los,reference_cases[c].sites,0,&out,&work)==TINYGAUSS_OK);CHECK(out.admissible_roots>0);double er=0,ev=0;for(size_t i=0;i<3;++i){double dr=out.position_m[i]-reference_cases[c].expected_position[i],dv=out.velocity_m_s[i]-reference_cases[c].expected_velocity[i];er+=dr*dr;ev+=dv*dv;CHECK(out.ranges_m[i]>0);}er=sqrt(er);ev=sqrt(ev);CHECK(er<1e-4&&ev<2e-5);max_r=fmax(max_r,er);max_v=fmax(max_v,ev);tinygauss_result selected;CHECK(tinygauss_estimate(reference_cases[c].mu,reference_cases[c].t1,reference_cases[c].t3,reference_cases[c].los,reference_cases[c].sites,out.selected_radius_m,&selected,&work)==TINYGAUSS_OK);CHECK(fabs(selected.selected_radius_m-out.selected_radius_m)<1e-6);}
 {tinygauss_result before=out;double parallel[]={1,0,0,1,0,0,1,0,0};CHECK(tinygauss_estimate(reference_cases[0].mu,-20,20,parallel,reference_cases[0].sites,0,&out,&work)==TINYGAUSS_GEOMETRY);CHECK(!memcmp(&out,&before,sizeof out));CHECK(tinygauss_estimate(0,-20,20,parallel,reference_cases[0].sites,0,&out,&work)==TINYGAUSS_INVALID_ARGUMENT);CHECK(tinygauss_estimate(reference_cases[0].mu,0,20,reference_cases[0].los,reference_cases[0].sites,0,&out,&work)==TINYGAUSS_INVALID_ARGUMENT);CHECK(tinygauss_estimate(reference_cases[0].mu,-20,20,NULL,reference_cases[0].sites,0,&out,&work)==TINYGAUSS_INVALID_ARGUMENT);}
 printf("tinygauss: 48 Vallado Gaussian/Gibbs oracle cases; max %.3g m / %.3g m/s\n",max_r,max_v);return 0;}
