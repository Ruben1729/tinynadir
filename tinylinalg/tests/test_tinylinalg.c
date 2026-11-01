/* SPDX-License-Identifier: MIT */
#include "tinylinalg/tinylinalg.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(c) do {if(!(c)){fprintf(stderr,"FAIL %s:%d: %s\n",__FILE__,__LINE__,#c);return 1;}}while(0)
#include "data/reference.h"
static int near(double a,double b,double tol){return fabs(a-b)<=tol*fmax(1,fabs(b));}
int main(void){
    double work[768],x[256],v[256],w[16],lu[256],lower[256];size_t pivots[16];
    double max_solve=0,max_lower=0,max_eigen=0;
    for(size_t c=0;c<sizeof reference_cases/sizeof reference_cases[0];++c){
        size_t n=reference_cases[c].n;const double *a=reference_cases[c].a,*b=reference_cases[c].b;
        CHECK(tinylinalg_solve(n,a,b,x,work,768)==TINYLINALG_OK);
        for(size_t i=0;i<n;++i){CHECK(near(x[i],reference_cases[c].x[i],2e-12));max_solve=fmax(max_solve,fabs(x[i]-reference_cases[c].x[i]));}
        CHECK(tinylinalg_lu(n,a,lu,pivots,work,768)==TINYLINALG_OK);
        CHECK(tinylinalg_lu_solve(n,lu,pivots,b,x,work,768)==TINYLINALG_OK);
        for(size_t i=0;i<n;++i)CHECK(near(x[i],reference_cases[c].x[i],2e-12));
        CHECK(tinylinalg_cholesky(n,a,lower,work,768)==TINYLINALG_OK);
        for(size_t i=0;i<n*n;++i){CHECK(near(lower[i],reference_cases[c].lower[i],2e-12));max_lower=fmax(max_lower,fabs(lower[i]-reference_cases[c].lower[i]));}
        CHECK(tinylinalg_eigen_symmetric(n,a,w,v,work,768)==TINYLINALG_OK);
        for(size_t i=0;i<n;++i){CHECK(near(w[i],reference_cases[c].eigenvalues[i],3e-12));max_eigen=fmax(max_eigen,fabs(w[i]-reference_cases[c].eigenvalues[i]));}
        for(size_t k=0;k<n;++k)for(size_t i=0;i<n;++i){double sum=0;for(size_t j=0;j<n;++j)sum+=a[i*n+j]*v[j*n+k];CHECK(near(sum,w[k]*v[i*n+k],3e-12));}
        for(size_t k=0;k<n;++k)for(size_t l=0;l<n;++l){double dot=0;for(size_t j=0;j<n;++j)dot+=v[j*n+k]*v[j*n+l];CHECK(near(dot,k==l?1:0,2e-12));}
        CHECK(tinylinalg_transpose(n,n,a,x,work,768)==TINYLINALG_OK);
        CHECK(tinylinalg_multiply(n,n,n,a,v,x,work,768)==TINYLINALG_OK);
        for(size_t i=0;i<n;++i)for(size_t j=0;j<n;++j)CHECK(near(x[i*n+j],v[i*n+j]*w[j],3e-12));
    }
    {const double a[]={0,2,1,3},b[]={4,7};CHECK(tinylinalg_solve(2,a,b,x,work,768)==TINYLINALG_OK);CHECK(near(x[0],1,1e-15)&&near(x[1],2,1e-15));}
    {const double a[]={1,2,2,4},b[]={1,2};x[0]=123;x[1]=456;CHECK(tinylinalg_solve(2,a,b,x,work,768)==TINYLINALG_SINGULAR);CHECK(x[0]==123&&x[1]==456);CHECK(tinylinalg_cholesky(2,a,lower,work,768)==TINYLINALG_NOT_POSITIVE_DEFINITE);}
    {const double a[]={1,2,0,4};CHECK(tinylinalg_eigen_symmetric(2,a,w,v,work,768)==TINYLINALG_INVALID_ARGUMENT);}
    {double a[]={NAN};CHECK(tinylinalg_solve(1,a,a,x,work,768)==TINYLINALG_INVALID_ARGUMENT);CHECK(tinylinalg_solve(0,a,a,x,work,768)==TINYLINALG_INVALID_ARGUMENT);CHECK(tinylinalg_solve(17,a,a,x,work,768)==TINYLINALG_INVALID_ARGUMENT);CHECK(tinylinalg_solve(1,a,a,x,NULL,0)==TINYLINALG_INVALID_ARGUMENT);}
    {const double a[]={DBL_MAX,DBL_MAX,0},b[]={0,0,0};CHECK(tinylinalg_unit3(a,x)==TINYLINALG_OK);CHECK(near(x[0],1/sqrt(2),1e-15));CHECK(tinylinalg_unit3(b,x)==TINYLINALG_SINGULAR);CHECK(tinylinalg_cross3(a,a,x)==TINYLINALG_ARITHMETIC_ERROR);}
    {double a[]={1,2,3},b[]={4,5,6};CHECK(tinylinalg_dot3(a,b)==32);CHECK(tinylinalg_cross3(a,b,a)==TINYLINALG_OK);CHECK(a[0]==-3&&a[1]==6&&a[2]==-3);}
    {double a[]={1,2,3,4};CHECK(tinylinalg_transpose(2,2,a,a,work,768)==TINYLINALG_OK);CHECK(a[1]==3&&a[2]==2);}
    printf("tinylinalg: 24 NumPy/SciPy cases; max solve %.3g, Cholesky %.3g, eigenvalue %.3g\n",max_solve,max_lower,max_eigen);
    return 0;
}
