/* SPDX-License-Identifier: MIT */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL %s:%d %s\n",__FILE__,__LINE__,#c); return 1; } } while (0)

#include "tinypolynomial/tinypolynomial.h"
#include "data/reference.h"
int main(void){tinypolynomial_complex roots[8],work[16];double maximum=0;
    for(size_t i=0;i<sizeof cases/sizeof cases[0];++i){const reference_case *c=&cases[i];int used[8]={0};CHECK(tinypolynomial_roots(c->degree,c->coefficients,roots,work,16)==0);for(unsigned j=0;j<c->degree;++j){size_t best=0;double distance=INFINITY;for(unsigned k=0;k<c->degree;++k)if(!used[k]){double d=hypot(roots[j].real-c->expected[k].real,roots[j].imaginary-c->expected[k].imaginary);if(d<distance){distance=d;best=k;}}used[best]=1;maximum=fmax(maximum,distance/fmax(1,hypot(c->expected[best].real,c->expected[best].imaginary)));CHECK(distance<2e-8*fmax(1,hypot(c->expected[best].real,c->expected[best].imaginary)));}}
    {double a[4]={1,0,0,0};CHECK(tinypolynomial_roots(3,a,roots,work,16)==0&&roots[0].real==0);a[0]=0;roots[0].real=123;CHECK(tinypolynomial_roots(3,a,roots,work,16)==TINYPOLYNOMIAL_INVALID&&roots[0].real==123);a[0]=1;a[1]=NAN;CHECK(tinypolynomial_roots(3,a,roots,work,16)==TINYPOLYNOMIAL_INVALID);CHECK(tinypolynomial_roots(9,a,roots,work,16)==TINYPOLYNOMIAL_INVALID);}
    printf("tinypolynomial: NumPy companion-root max scaled difference %.12g\n",maximum);return 0;}
