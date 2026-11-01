/* SPDX-License-Identifier: MIT */
#include "tinylambert/tinylambert.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "data/reference.h"
#define CHECK(x) do { if (!(x)) {fprintf(stderr,"line %d: %s\n",__LINE__,#x);return 1;} } while(0)
int main(void)
{
    size_t i,j;double maximum=0;
    for(i=0;i<sizeof reference/sizeof reference[0];++i){const double *r=reference[i];tinylambert_solution out;double e1=0,e2=0,r1=hypot(hypot(r[1],r[2]),r[3]),r2=hypot(hypot(r[4],r[5]),r[6]);
        CHECK(tinylambert_solve(r[0],r+1,r+4,r[7],(tinylambert_path)(int)r[8],&out)==TINYLAMBERT_OK);
        for(j=0;j<3;++j){double a=fabs(out.departure_m_s[j]-r[9+j]),b=fabs(out.arrival_m_s[j]-r[12+j]);if(a>maximum)maximum=a;if(b>maximum)maximum=b;if(a>=2e-6||b>=2e-6)fprintf(stderr,"case %zu: %.9g %.9g m/s\n",i,a,b);CHECK(a<2e-6);CHECK(b<2e-6);e1+=.5*out.departure_m_s[j]*out.departure_m_s[j];e2+=.5*out.arrival_m_s[j]*out.arrival_m_s[j];}
        e1-=r[0]/r1;e2-=r[0]/r2;CHECK(fabs(e1-e2)<2e-6);CHECK(out.iterations<=160);
    }
    {
        double a[3]={7000000,0,0},b[3]={0,7000000,0};tinylambert_solution out={{1,2,3},{4,5,6},7},saved=out;
        CHECK(tinylambert_solve(3.986e14,a,a,1000,TINYLAMBERT_SHORT_WAY,&out)==TINYLAMBERT_GEOMETRY);CHECK(!memcmp(&out,&saved,sizeof out));
        CHECK(tinylambert_solve(3.986e14,a,b,-1,TINYLAMBERT_SHORT_WAY,&out)==TINYLAMBERT_INVALID);CHECK(!memcmp(&out,&saved,sizeof out));
        CHECK(tinylambert_solve(3.986e14,a,b,1e-10,TINYLAMBERT_SHORT_WAY,&out)==TINYLAMBERT_RANGE);
    }
    printf("Lambert: 18 Vallado short/long zero-revolution cases; max component error %.9g m/s\n",maximum);return 0;
}
