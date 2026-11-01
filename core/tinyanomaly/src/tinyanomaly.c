/* SPDX-License-Identifier: MIT */
#include "tinyanomaly/tinyanomaly.h"
#include <math.h>
#include <float.h>
#define TA_PI 3.14159265358979323846264338327950288

/* Preserve the small cubic term near the parabolic limit. */
static double difference(double x, int hyperbolic)
{
    if (fabs(x)<0.5) {
        double term=x*x*x/6.0, sum=term;
        unsigned k;
        for (k=2; k<=12; ++k) {
            term *= (hyperbolic ? 1.0 : -1.0)*x*x/
                    ((double)(2*k)*(double)(2*k+1));
            sum += term;
        }
        return sum;
    }
    return hyperbolic ? sinh(x)-x : x-sin(x);
}

tinyanomaly_status tinyanomaly_from_mean(double e, double mean, tinyanomaly_solution *out)
{
    tinyanomaly_solution result;
    double lo, hi, x, m;
    unsigned iteration;
    if (!out || !isfinite(e) || e<0.0 || e>1000.0 ||
        !isfinite(mean) || fabs(mean)>1e12) return TINYANOMALY_INVALID;
    if (e==1.0) {
        result.anomaly_rad=2.0*sinh(asinh(1.5*mean)/3.0);
        result.true_anomaly_rad=2.0*atan(result.anomaly_rad);
        *out=result;
        return TINYANOMALY_OK;
    }
    m=e<1.0 ? remainder(mean,2.0*TA_PI) : mean;
    if (m==0.0) { result.anomaly_rad=0.0;result.true_anomaly_rad=0.0;*out=result;return TINYANOMALY_OK; }
    lo=e<1.0 ? -TA_PI : -(log1p(2.0*fabs(m)/e)+2.0);
    hi=-lo;
    x=e<1.0 ? m : asinh(m/e);
    for (iteration=0; iteration<100; ++iteration) {
        double f,derivative,next;
        if (e<1.0) {
            double s=sin(0.5*x);
            f=(1.0-e)*x+e*difference(x,0)-m;
            derivative=(1.0-e)+2.0*e*s*s;
        } else {
            double s=sinh(0.5*x);
            f=(e-1.0)*x+e*difference(x,1)-m;
            derivative=(e-1.0)+2.0*e*s*s;
        }
        if (f==0.0) break;
        if (f>0.0) hi=x; else lo=x;
        next=x-f/derivative;
        if (!(next>lo && next<hi)) next=lo+0.5*(hi-lo);
        if (next==x || hi-lo<=4.0*DBL_EPSILON*(fabs(x)+1e-6)) break;
        x=next;
    }
    if (iteration==100) return TINYANOMALY_NO_CONVERGENCE;
    result.anomaly_rad=x;
    if (e<1.0)
        result.true_anomaly_rad=2.0*atan2(sqrt(1.0+e)*sin(0.5*x),
                                         sqrt(1.0-e)*cos(0.5*x));
    else result.true_anomaly_rad=2.0*atan(sqrt((e+1.0)/(e-1.0))*tanh(0.5*x));
    if (!isfinite(result.anomaly_rad) || !isfinite(result.true_anomaly_rad))
        return TINYANOMALY_NO_CONVERGENCE;
    *out=result;
    return TINYANOMALY_OK;
}

tinyanomaly_status tinyanomaly_from_true(double e, double nu, double *mean, double *anomaly)
{
    double a,m;
    if (!mean || !anomaly || !isfinite(e) || e<0.0 || e>1000.0 ||
        !isfinite(nu) || fabs(nu)>TA_PI) return TINYANOMALY_INVALID;
    if (e<1.0) {
        a=2.0*atan2(sqrt(1.0-e)*sin(nu/2.0),sqrt(1.0+e)*cos(nu/2.0));
        m=(1.0-e)*a+e*difference(a,0);
    } else if (e==1.0) {
        if (fabs(nu)>=TA_PI) return TINYANOMALY_INVALID;
        a=tan(nu/2.0);m=a+a*a*a/3.0;
    } else {
        double arg;
        if (fabs(nu)>=acos(-1.0/e)) return TINYANOMALY_INVALID;
        arg=sqrt((e-1.0)/(e+1.0))*tan(nu/2.0);
        a=2.0*atanh(arg);m=(e-1.0)*a+e*difference(a,1);
    }
    if (!isfinite(a) || !isfinite(m) || fabs(m)>1e12) return TINYANOMALY_INVALID;
    *mean=m;*anomaly=a;
    return TINYANOMALY_OK;
}
