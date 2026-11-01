#include "tinymoon/tinymoon.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(x) do{if(!(x)){fprintf(stderr,"%s:%d %s\n",__FILE__,__LINE__,#x);exit(1);}}while(0)
#include "data/reference.h"
int main(void) {
    tinymoon_state state,saved;size_t i,k;double maxr=0,maxv=0,direction[3],position[3]={7000000,1000000,-2000000};
    for(i=0;i<sizeof references/sizeof references[0];i++) {
        CHECK(tinymoon_geocentric(references[i][0],references[i][1],&state)==TINYMOON_OK);
        for(k=0;k<3;k++) {double r=fabs(state.position_m[k]-references[i][2+k]),v=fabs(state.velocity_m_s[k]-references[i][5+k]);CHECK(r<.001);CHECK(v<1e-8);if(r>maxr)maxr=r;if(v>maxv)maxv=v;}
        CHECK(tinymoon_direction(&state,position,direction)==TINYMOON_OK);CHECK(fabs(hypot(hypot(direction[0],direction[1]),direction[2])-1)<3e-16);
    }
    saved=state;
    CHECK(tinymoon_geocentric(NAN,0,&state)==TINYMOON_INVALID);CHECK(!memcmp(&state,&saved,sizeof state));
    CHECK(tinymoon_geocentric(2451545,36526,&state)==TINYMOON_RANGE);CHECK(!memcmp(&state,&saved,sizeof state));
    CHECK(tinymoon_geocentric(2451545,0,NULL)==TINYMOON_INVALID);
    direction[0]=1;direction[1]=2;direction[2]=3;
    CHECK(tinymoon_direction(&state,state.position_m,direction)==TINYMOON_GEOMETRY);CHECK(direction[0]==1&&direction[1]==2&&direction[2]==3);
    CHECK(tinymoon_direction(NULL,position,direction)==TINYMOON_INVALID);
    printf("tinymoon: %zu ERFA reference states and direction/failure cases pass; max %.9g m / %.9g m/s\n",sizeof references/sizeof references[0],maxr,maxv);
    return 0;
}
