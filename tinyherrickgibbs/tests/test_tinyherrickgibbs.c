/* SPDX-License-Identifier: MIT */
#include "tinyherrickgibbs/tinyherrickgibbs.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "data/reference.h"
#define CHECK(x) do { if (!(x)) {fprintf(stderr,"line %d: %s\n",__LINE__,#x);return 1;} } while(0)
int main(void)
{
    size_t i,j;double maximum=0;
    for(i=0;i<sizeof reference/sizeof reference[0];++i) {
        tinyherrickgibbs_result out;
        CHECK(tinyherrickgibbs_estimate(reference[i][0],reference[i]+3,reference[i]+6,reference[i]+9,reference[i][1],reference[i][2],.05,.01,&out)==TINYHERRICKGIBBS_OK);
        for(j=0;j<3;++j){double error=fabs(out.velocity_m_s[j]-reference[i][12+j]);if(error>maximum)maximum=error;CHECK(error<2e-8);}
        CHECK(fabs(out.angle12_rad-reference[i][15])<2e-10);CHECK(fabs(out.angle23_rad-reference[i][16])<2e-10);CHECK(fabs(out.coplanarity_rad-reference[i][17])<2e-12);
    }
    {
        double r[3]={1,0,0};tinyherrickgibbs_result out={{1,2,3},4,5,6},saved=out;
        CHECK(tinyherrickgibbs_estimate(1,r,r,r,1,1,.1,.1,&out)==TINYHERRICKGIBBS_GEOMETRY);CHECK(!memcmp(&out,&saved,sizeof out));
        CHECK(tinyherrickgibbs_estimate(1,r,r,r,0,1,.1,.1,&out)==TINYHERRICKGIBBS_INVALID);CHECK(!memcmp(&out,&saved,sizeof out));
        CHECK(tinyherrickgibbs_estimate(1,r,r,r,1,NAN,.1,.1,&out)==TINYHERRICKGIBBS_INVALID);
    }
    printf("Herrick-Gibbs: 15 frozen Vallado cases; max component error %.9g m/s\n",maximum);
    return 0;
}
