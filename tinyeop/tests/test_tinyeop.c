#include "tinyeop/tinyeop.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(x) do{if(!(x)){fprintf(stderr,"%s:%d %s\n",__FILE__,__LINE__,#x);exit(1);}}while(0)
#include "data/reference.h"
int main(void) {
    size_t i,k;double maximum=0,v[6],xys[3];tinyeop_result out,before;tinyeop_record bad[8];
    tinyeop_xys_record points[8];
    for(i=0;i<8;i++){points[i].mjd_tt=input_records[i].mjd_utc;points[i].x_rad=input_records[i].value[0];points[i].y_rad=input_records[i].value[4];points[i].s_rad=input_records[i].value[5];}
    for(i=0;i<sizeof references/sizeof references[0];i++) {
        CHECK(tinyeop_sample(input_records,8,references[i].t,37,(tinyeop_method)references[i].method,&out)==TINYEOP_OK);
        v[0]=out.xp_rad;v[1]=out.yp_rad;v[2]=out.dut1_s;v[3]=out.lod_s;v[4]=out.dx_rad;v[5]=out.dy_rad;
        for(k=0;k<6;k++){double e=fabs(v[k]-references[i].v[k]);CHECK(e<2e-13);if(e>maximum)maximum=e;}
        CHECK(tinyeop_sample_xys(points,8,references[i].t,(tinyeop_method)references[i].method,xys)==TINYEOP_OK);
        CHECK(fabs(xys[0]-references[i].v[0])<1e-19);CHECK(fabs(xys[1]-references[i].v[4])<1e-21);CHECK(fabs(xys[2]-references[i].v[5])<1e-21);
    }
    CHECK(tinyeop_sample(input_records,8,60003.5,37,TINYEOP_LINEAR,&out)==TINYEOP_OK);before=out;
    CHECK(tinyeop_sample(input_records,8,60003.5,38,TINYEOP_LINEAR,&out)==TINYEOP_OK);CHECK(fabs(out.dut1_s-before.dut1_s-1)<1e-14);
    before=out;
    CHECK(tinyeop_sample(input_records,8,59999,37,TINYEOP_LINEAR,&out)==TINYEOP_RANGE);CHECK(!memcmp(&out,&before,sizeof out));
    CHECK(tinyeop_sample(input_records,8,60008,37,TINYEOP_LINEAR,&out)==TINYEOP_RANGE);
    CHECK(tinyeop_sample(input_records,3,60001,37,TINYEOP_CUBIC,&out)==TINYEOP_RANGE);
    CHECK(tinyeop_sample(input_records,0,60001,37,TINYEOP_LINEAR,&out)==TINYEOP_RANGE);
    CHECK(tinyeop_sample(NULL,8,60001,37,TINYEOP_LINEAR,&out)==TINYEOP_INVALID);
    CHECK(tinyeop_sample(input_records,8,NAN,37,TINYEOP_LINEAR,&out)==TINYEOP_INVALID);
    memcpy(bad,input_records,sizeof bad);bad[2].mjd_utc=bad[1].mjd_utc;
    CHECK(tinyeop_sample(bad,8,60001,37,TINYEOP_LINEAR,&out)==TINYEOP_TABLE);
    memcpy(bad,input_records,sizeof bad);bad[3].value[0]=NAN;
    CHECK(tinyeop_sample(bad,8,60001,37,TINYEOP_LINEAR,&out)==TINYEOP_TABLE);CHECK(!memcmp(&out,&before,sizeof out));
    printf("tinyeop: %zu NumPy comparisons, XYS, leap-offset discontinuity and failure checks pass; max error %.9g\n",sizeof references/sizeof references[0],maximum);
    return 0;
}
