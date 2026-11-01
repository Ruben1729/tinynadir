/* SPDX-License-Identifier: MIT */
#include "tinystumpff/tinystumpff.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "data/reference.h"
#define CHECK(x) do { if (!(x)) { fprintf(stderr,"line %d: %s\n",__LINE__,#x);return 1;} } while(0)
int main(void)
{
    size_t i,k;
    double maximum=0.0;
    for(i=0;i<sizeof reference/sizeof reference[0];++i) {
        tinystumpff_values out;
        CHECK(tinystumpff_evaluate(reference[i][0],&out)==TINYSTUMPFF_OK);
        for(k=0;k<4;++k) {
            double scale=fmax(1.0,fabs(reference[i][k+1]));
            double error=fabs(out.c[k]-reference[i][k+1])/scale;
            if(error>maximum)maximum=error;
            CHECK(error<2e-13);
        }
    }
    {
        tinystumpff_values out={{11,12,13,14}},saved=out;
        CHECK(tinystumpff_evaluate(NAN,&out)==TINYSTUMPFF_INVALID);
        CHECK(tinystumpff_evaluate(-500000.,&out)==TINYSTUMPFF_RANGE);
        CHECK(!memcmp(&out,&saved,sizeof out));
        CHECK(tinystumpff_evaluate(0.,NULL)==TINYSTUMPFF_INVALID);
    }
    printf("Stumpff: %zu reference cases; max scaled error %.9g\n",sizeof reference/sizeof reference[0],maximum);
    return 0;
}
