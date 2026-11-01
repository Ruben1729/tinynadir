/* SPDX-License-Identifier: MIT */
#include "tinygibbs/tinygibbs.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "data/reference.h"
#define CHECK(x) do { if (!(x)) {fprintf(stderr,"line %d: %s\n",__LINE__,#x);return 1;} } while(0)
int main(void)
{
    size_t i,j;double maximum=0;
    for(i=0;i<sizeof reference/sizeof reference[0];++i) {
        tinygibbs_result out;
        CHECK(tinygibbs_estimate(reference[i][0],reference[i]+1,reference[i]+4,reference[i]+7,.01,&out)==TINYGIBBS_OK);
        for(j=0;j<3;++j){double error=fabs(out.velocity_m_s[j]-reference[i][10+j]);if(error>maximum)maximum=error;CHECK(error<2e-6);}
        CHECK(fabs(out.angle12_rad-reference[i][13])<2e-12);CHECK(fabs(out.angle23_rad-reference[i][14])<2e-12);CHECK(fabs(out.coplanarity_rad-reference[i][15])<2e-12);
    }
    {
        double r1[3]={1,0,0},r2[3]={2,0,0},r3[3]={3,0,0};tinygibbs_result out={{1,2,3},4,5,6},saved=out;
        CHECK(tinygibbs_estimate(1,r1,r2,r3,.01,&out)==TINYGIBBS_GEOMETRY);CHECK(!memcmp(&out,&saved,sizeof out));
        r2[0]=0;r2[1]=1;r3[0]=0;r3[2]=1;
        CHECK(tinygibbs_estimate(1,r1,r2,r3,.01,&out)==TINYGIBBS_GEOMETRY);CHECK(!memcmp(&out,&saved,sizeof out));
        CHECK(tinygibbs_estimate(-1,r1,r2,r3,.01,&out)==TINYGIBBS_INVALID);
    }
    printf("Gibbs: 15 frozen Vallado cases; max component error %.9g m/s\n",maximum);
    return 0;
}
