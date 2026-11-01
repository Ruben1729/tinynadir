/* SPDX-License-Identifier: MIT */
#include "tinypolynomial/tinypolynomial.h"
#include <complex.h>
#include <float.h>
#include <math.h>
#include <string.h>
#define TP_PI 3.14159265358979323846264338327950288
static double complex get(tinypolynomial_complex a){return a.real+I*a.imaginary;}
static tinypolynomial_complex put(double complex a){tinypolynomial_complex r={creal(a),cimag(a)};return r;}
static int finite_complex(double complex a){return isfinite(creal(a))&&isfinite(cimag(a));}
tinypolynomial_status tinypolynomial_roots(unsigned n,const double *a,tinypolynomial_complex *out,tinypolynomial_complex *w,size_t count){
    double b[9],logscale=-INFINITY,scale_value;tinypolynomial_complex result[8];int all_zero=1,converged=0;
    if(!n||n>8||!a||!out||!w||count<2*n||!isfinite(a[0])||a[0]==0)return TINYPOLYNOMIAL_INVALID;
    for(unsigned i=1;i<=n;++i){if(!isfinite(a[i]))return TINYPOLYNOMIAL_INVALID;if(a[i]!=0){all_zero=0;logscale=fmax(logscale,(log(fabs(a[i]))-log(fabs(a[0])))/i);}}
    if(all_zero){for(unsigned i=0;i<n;++i)result[i]=put(0);memcpy(out,result,n*sizeof *out);return TINYPOLYNOMIAL_OK;}
    if(!isfinite(logscale)||logscale>log(DBL_MAX)||logscale<log(DBL_MIN))return TINYPOLYNOMIAL_ARITHMETIC;
    scale_value=exp(logscale);b[0]=1;
    for(unsigned i=1;i<=n;++i)b[i]=a[i]==0?0:copysign(1,a[i])*copysign(1,a[0])*exp(log(fabs(a[i]))-log(fabs(a[0]))-i*logscale);
    if(n==1){w[0]=put(-b[1]);converged=1;}
    else if(n==2){double complex root=csqrt(b[1]*b[1]-4*b[2]+0.0*I),first=-.5*(b[1]+(b[1]>=0?root:-root));w[0]=put(first);w[1]=put(first==0?-b[1]:b[2]/first);converged=1;}
    else {
        double radius=1;for(unsigned i=1;i<=n;++i)radius=fmax(radius,1+fabs(b[i]));
        for(unsigned i=0;i<n;++i){double phase=2*TP_PI*((double)i+.371)/(double)n;w[i]=put(radius*(cos(phase)+I*sin(phase)));}
        for(unsigned iteration=0;iteration<256;++iteration){double largest=0;
            for(unsigned i=0;i<n;++i){double complex z=get(w[i]),p=1,dp=0,sum=0,correction,denominator;
                for(unsigned j=1;j<=n;++j){dp=dp*z+p;p=p*z+b[j];}
                for(unsigned j=0;j<n;++j)if(i!=j){double complex difference=z-get(w[j]);if(cabs(difference)<DBL_MIN)return TINYPOLYNOMIAL_NO_CONVERGENCE;sum+=1/difference;}
                denominator=dp-p*sum;if(cabs(denominator)<DBL_MIN)return TINYPOLYNOMIAL_NO_CONVERGENCE;
                correction=p/denominator;if(!finite_complex(correction))return TINYPOLYNOMIAL_NO_CONVERGENCE;
                if(cabs(correction)>1+cabs(z))correction*= (1+cabs(z))/cabs(correction);
                w[n+i]=put(z-correction);largest=fmax(largest,cabs(correction)/(1+cabs(z)));
            }
            memcpy(w,w+n,n*sizeof *w);
            if(largest<4e-14){converged=1;break;}
        }
    }
    if(!converged)return TINYPOLYNOMIAL_NO_CONVERGENCE;
    for(unsigned i=0;i<n;++i){double complex z=get(w[i]),p=1;double bound=1;for(unsigned j=1;j<=n;++j){p=p*z+b[j];bound=bound*cabs(z)+fabs(b[j]);}if(cabs(p)>2e-11*fmax(DBL_MIN,bound))return TINYPOLYNOMIAL_NO_CONVERGENCE;z*=scale_value;if(!finite_complex(z))return TINYPOLYNOMIAL_ARITHMETIC;result[i]=put(z);}
    for(unsigned i=0;i<n;++i)for(unsigned j=i+1;j<n;++j)if(result[j].real<result[i].real||(result[j].real==result[i].real&&result[j].imaginary<result[i].imaginary)){tinypolynomial_complex temp=result[i];result[i]=result[j];result[j]=temp;}
    memcpy(out,result,n*sizeof *out);return TINYPOLYNOMIAL_OK;
}
