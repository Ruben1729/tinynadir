/* SPDX-License-Identifier: MIT */
#include "tinystumpff/tinystumpff.h"
#include <math.h>
#include <float.h>

tinystumpff_status tinystumpff_evaluate(double z, tinystumpff_values *output)
{
    tinystumpff_values result;
    unsigned k;
    if (!output || !isfinite(z)) return TINYSTUMPFF_INVALID;
    if (fabs(z) < 0.5) {
        const double initial[4] = {1.0, 1.0, 0.5, 1.0/6.0};
        for (k=0; k<4; ++k) {
            double term=initial[k], sum=term;
            unsigned n;
            for (n=0; n<40; ++n) {
                double denominator=(double)(2*n+k+1)*(double)(2*n+k+2);
                double next;
                term *= -z/denominator;
                next=sum+term;
                if (next==sum) break;
                sum=next;
            }
            result.c[k]=sum;
        }
    } else if (z>0.0) {
        double u=sqrt(z), half_sine=sin(0.5*u);
        result.c[0]=cos(u);
        result.c[1]=sin(u)/u;
        result.c[2]=2.0*half_sine*half_sine/z;
        result.c[3]=(1.0-result.c[1])/z;
    } else {
        double u=sqrt(-z), half_sine;
        if (u>700.0) return TINYSTUMPFF_RANGE;
        half_sine=sinh(0.5*u);
        result.c[0]=cosh(u);
        result.c[1]=sinh(u)/u;
        result.c[2]=2.0*half_sine*half_sine/(-z);
        result.c[3]=(result.c[1]-1.0)/(-z);
    }
    for (k=0; k<4; ++k)
        if (!isfinite(result.c[k])) return TINYSTUMPFF_RANGE;
    *output=result;
    return TINYSTUMPFF_OK;
}
