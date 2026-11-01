#include "tinydtm/tinydtm.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(x) do{if(!(x)){fprintf(stderr,"%s:%d %s\n",__FILE__,__LINE__,#x);exit(1);}}while(0)
#include "data/reference.h"
int main(void) {
    tinydtm_workspace workspace,other;tinydtm_output out,saved,result;tinydtm_input in;size_t i,k;double max_relative=0;
    for(i=0;i<sizeof references/sizeof references[0];i++) {
        double values[9],sum=0;
        CHECK(tinydtm_evaluate(&references[i].in,&workspace,&out)==TINYDTM_OK);
        values[0]=out.density_kg_m3;values[1]=out.temperature_k;values[2]=out.exospheric_temperature_k;
        for(k=0;k<6;k++){values[3+k]=out.number_density_m3[k];sum+=out.partial_mass_density_kg_m3[k];}
        for(k=0;k<9;k++){double e=fabs(values[k]-references[i].expected[k]);CHECK(e<1e-24+2e-11*fabs(references[i].expected[k]));if(references[i].expected[k]!=0&&e/fabs(references[i].expected[k])>max_relative)max_relative=e/fabs(references[i].expected[k]);}
        CHECK(fabs(sum-out.density_kg_m3)<1e-24+1e-13*out.density_kg_m3);
        CHECK(tinydtm_evaluate(&references[i].in,&other,&result)==TINYDTM_OK);CHECK(!memcmp(&out,&result,sizeof out));
    }
    saved=out;in=references[0].in;in.altitude_m=119999;
    CHECK(tinydtm_evaluate(&in,&workspace,&out)==TINYDTM_RANGE);CHECK(!memcmp(&out,&saved,sizeof out));
    in=references[0].in;in.f107=NAN;CHECK(tinydtm_evaluate(&in,&workspace,&out)==TINYDTM_INVALID);CHECK(!memcmp(&out,&saved,sizeof out));
    in=references[0].in;in.kp_3h=10;CHECK(tinydtm_evaluate(&in,&workspace,&out)==TINYDTM_RANGE);
    CHECK(tinydtm_evaluate(NULL,&workspace,&out)==TINYDTM_INVALID);CHECK(tinydtm_evaluate(&references[0].in,NULL,&out)==TINYDTM_INVALID);
    printf("tinydtm: %zu Orekit DTM2000 comparisons pass; max relative error %.9g; workspace %zu bytes\n",sizeof references/sizeof references[0],max_relative,sizeof workspace);
    return 0;
}
