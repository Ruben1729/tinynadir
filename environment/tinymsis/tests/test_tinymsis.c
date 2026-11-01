#include "tinymsis/tinymsis.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(x) do{if(!(x)){fprintf(stderr,"%s:%d %s\n",__FILE__,__LINE__,#x);exit(1);}}while(0)
#include "data/reference.h"
int main(void) {
    tinymsis_workspace w,other;tinymsis_output out,saved,r;size_t i,k;double maximum=0;tinymsis_input bad;int switches[24];
    for(i=0;i<sizeof references/sizeof references[0];i++) {
        double value[11];
        CHECK(tinymsis_evaluate(&references[i].input,references[i].switches,&w,&out)==TINYMSIS_OK);
        for(k=0;k<8;k++)value[k]=out.number_density_m3[k];value[8]=out.density_kg_m3;value[9]=out.exospheric_temperature_k;value[10]=out.temperature_k;
        for(k=0;k<11;k++){double e=fabs(value[k]-references[i].expected[k]),floor=k==8?1e-30:1e-8;CHECK(e<=floor+1e-10*fabs(references[i].expected[k]));if(references[i].expected[k]&&e/fabs(references[i].expected[k])>maximum)maximum=e/fabs(references[i].expected[k]);}
        CHECK(tinymsis_evaluate(&references[i].input,references[i].switches,&other,&r)==TINYMSIS_OK);CHECK(!memcmp(&out,&r,sizeof out));
    }
    CHECK(tinymsis_evaluate(&references[0].input,NULL,&w,&out)==TINYMSIS_OK);saved=out;
    tinymsis_default_switches(switches);for(i=0;i<24;i++)CHECK(switches[i]==1);
    CHECK(tinymsis_evaluate(&references[1].input,NULL,&w,&r)==TINYMSIS_OK);
    CHECK(tinymsis_evaluate(&references[0].input,NULL,&w,&out)==TINYMSIS_OK);CHECK(!memcmp(&out,&saved,sizeof out));
    bad=references[0].input;bad.altitude_m=-1;CHECK(tinymsis_evaluate(&bad,NULL,&w,&out)==TINYMSIS_RANGE);CHECK(!memcmp(&out,&saved,sizeof out));
    bad=references[0].input;bad.f107=NAN;CHECK(tinymsis_evaluate(&bad,NULL,&w,&out)==TINYMSIS_INVALID);CHECK(!memcmp(&out,&saved,sizeof out));
    bad=references[0].input;bad.ap[4]=-1;CHECK(tinymsis_evaluate(&bad,NULL,&w,&out)==TINYMSIS_RANGE);
    switches[0]=0;CHECK(tinymsis_evaluate(&references[0].input,switches,&w,&out)==TINYMSIS_RANGE);CHECK(!memcmp(&out,&saved,sizeof out));
    CHECK(tinymsis_evaluate(NULL,NULL,&w,&out)==TINYMSIS_INVALID);CHECK(tinymsis_evaluate(&references[0].input,NULL,NULL,&out)==TINYMSIS_INVALID);
    printf("tinymsis: %zu original-model comparisons pass; max relative error %.9g; caller workspace %zu bytes\n",sizeof references/sizeof references[0],maximum,sizeof w);
    return 0;
}
