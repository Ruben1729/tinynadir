#include "tinyframes/tinyframes.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(x) do{if(!(x)){fprintf(stderr,"%s:%d %s\n",__FILE__,__LINE__,#x);exit(1);}}while(0)
#include "data/reference.h"
int main(void) {
    size_t n;int i,j,k;double a[3][3],maximum=0,r[3],v[3],rr[3],vv[3],before[3][3];tinyframes_epoch bad;
    for(n=0;n<sizeof references/sizeof references[0];n++) {
        CHECK(tinyframes_matrix(&epochs[references[n].epoch],(tinyframes_frame)references[n].from,(tinyframes_frame)references[n].to,(tinyframes_model)references[n].model,a)==TINYFRAMES_OK);
        for(i=0;i<3;i++)for(j=0;j<3;j++){double e=fabs(a[i][j]-references[n].m[3*i+j]);CHECK(e<5e-14);if(e>maximum)maximum=e;}
        for(i=0;i<3;i++)for(j=0;j<3;j++){double s=0;for(k=0;k<3;k++)s+=a[i][k]*a[j][k];CHECK(fabs(s-(i==j))<5e-15);}
    }
    for(n=0;n<sizeof state_cases/sizeof state_cases[0];n++) {
        const tinyframes_epoch *epoch=&epochs[state_cases[n].epoch];
        tinyframes_frame from=(tinyframes_frame)state_cases[n].from,to=(tinyframes_frame)state_cases[n].to;tinyframes_model model=(tinyframes_model)state_cases[n].model;
        CHECK(tinyframes_transform(epoch,from,to,model,state_cases[n].r,state_cases[n].v,r,v)==TINYFRAMES_OK);
        for(i=0;i<3;i++){CHECK(fabs(r[i]-state_cases[n].rr[i])<1e-6);CHECK(fabs(v[i]-state_cases[n].vv[i])<1e-6);}
        CHECK(tinyframes_transform(epoch,to,from,model,r,v,rr,vv)==TINYFRAMES_OK);
        for(i=0;i<3;i++){CHECK(fabs(rr[i]-state_cases[n].r[i])<5e-8);CHECK(fabs(vv[i]-state_cases[n].v[i])<1e-6);}
    }
    bad=epochs[0];bad.xp_rad=NAN;memcpy(before,a,sizeof before);
    CHECK(tinyframes_matrix(&bad,TINYFRAMES_GCRF,TINYFRAMES_ITRF,TINYFRAMES_IAU2006A,a)==TINYFRAMES_INVALID);CHECK(!memcmp(before,a,sizeof a));
    bad=epochs[0];bad.tt[0]=100;
    CHECK(tinyframes_matrix(&bad,TINYFRAMES_GCRF,TINYFRAMES_ITRF,TINYFRAMES_IAU2006A,a)==TINYFRAMES_RANGE);CHECK(!memcmp(before,a,sizeof a));
    CHECK(tinyframes_matrix(NULL,TINYFRAMES_GCRF,TINYFRAMES_ITRF,TINYFRAMES_IAU2006A,a)==TINYFRAMES_INVALID);
    CHECK(tinyframes_matrix(&epochs[0],(tinyframes_frame)99,TINYFRAMES_ITRF,TINYFRAMES_IAU2006A,a)==TINYFRAMES_RANGE);
    CHECK(tinyframes_matrix(&epochs[0],(tinyframes_frame)-1,TINYFRAMES_ITRF,TINYFRAMES_IAU2006A,a)==TINYFRAMES_RANGE);CHECK(!memcmp(before,a,sizeof a));
    CHECK(tinyframes_matrix(&epochs[0],TINYFRAMES_GCRF,(tinyframes_frame)-1,TINYFRAMES_IAU2006A,a)==TINYFRAMES_RANGE);CHECK(!memcmp(before,a,sizeof a));
    CHECK(tinyframes_matrix(&epochs[0],TINYFRAMES_GCRF,TINYFRAMES_ITRF,(tinyframes_model)-1,a)==TINYFRAMES_RANGE);CHECK(!memcmp(before,a,sizeof a));
    printf("tinyframes: %zu ERFA matrices, %zu state cases and inverse/invalid checks pass; max matrix error %.9g\n",sizeof references/sizeof references[0],sizeof state_cases/sizeof state_cases[0],maximum);
    return 0;
}
