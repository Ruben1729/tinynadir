/* SPDX-License-Identifier: MIT */
#include "tinylinalg/tinylinalg.h"
#include <float.h>
#include <math.h>
#include <string.h>

static int valid_size(size_t n) { return n > 0 && n <= TINYLINALG_MAX_DIM; }
static int finite_array(const double *a, size_t n) {
    if (!a) return 0;
    for (size_t i = 0; i < n; ++i) if (!isfinite(a[i])) return 0;
    return 1;
}
double tinylinalg_dot3(const double a[3], const double b[3]) {
    if (!finite_array(a, 3) || !finite_array(b, 3)) return NAN;
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
double tinylinalg_norm3(const double v[3]) {
    if (!finite_array(v, 3)) return NAN;
    return hypot(hypot(v[0], v[1]), v[2]);
}
tinylinalg_status tinylinalg_cross3(const double a[3], const double b[3], double out[3]) {
    double r[3];
    if (!out || !finite_array(a, 3) || !finite_array(b, 3)) return TINYLINALG_INVALID_ARGUMENT;
    r[0]=a[1]*b[2]-a[2]*b[1];r[1]=a[2]*b[0]-a[0]*b[2];r[2]=a[0]*b[1]-a[1]*b[0];
    if (!finite_array(r, 3)) return TINYLINALG_ARITHMETIC_ERROR;
    memcpy(out,r,sizeof r);return TINYLINALG_OK;
}
tinylinalg_status tinylinalg_unit3(const double v[3], double out[3]) {
    double r[3], scale, norm;
    if (!out || !finite_array(v, 3)) return TINYLINALG_INVALID_ARGUMENT;
    scale=fmax(fabs(v[0]),fmax(fabs(v[1]),fabs(v[2])));
    if (scale==0) return TINYLINALG_SINGULAR;
    for(size_t i=0;i<3;++i)r[i]=v[i]/scale;
    norm=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
    for(size_t i=0;i<3;++i)r[i]/=norm;
    memcpy(out,r,sizeof r);return TINYLINALG_OK;
}
tinylinalg_status tinylinalg_multiply(size_t rows,size_t inner,size_t cols,
    const double *a,const double *b,double *out,double *work,size_t count) {
    if(!valid_size(rows)||!valid_size(inner)||!valid_size(cols)||!out||!work||count<rows*cols||
       !finite_array(a,rows*inner)||!finite_array(b,inner*cols))return TINYLINALG_INVALID_ARGUMENT;
    for(size_t i=0;i<rows;++i)for(size_t j=0;j<cols;++j){
        double x=0;for(size_t k=0;k<inner;++k)x+=a[i*inner+k]*b[k*cols+j];
        if(!isfinite(x))return TINYLINALG_ARITHMETIC_ERROR;
        work[i*cols+j]=x;
    }
    memcpy(out,work,rows*cols*sizeof *out);return TINYLINALG_OK;
}
tinylinalg_status tinylinalg_transpose(size_t rows,size_t cols,const double *a,
    double *out,double *work,size_t count) {
    if(!valid_size(rows)||!valid_size(cols)||!out||!work||count<rows*cols||!finite_array(a,rows*cols))return TINYLINALG_INVALID_ARGUMENT;
    for(size_t i=0;i<rows;++i)for(size_t j=0;j<cols;++j)work[j*rows+i]=a[i*cols+j];
    memcpy(out,work,rows*cols*sizeof *out);return TINYLINALG_OK;
}
tinylinalg_status tinylinalg_lu(size_t n,const double *a,double *lu,
    size_t *pivots,double *work,size_t count) {
    size_t p[TINYLINALG_MAX_DIM];
    if(!valid_size(n)||!lu||!pivots||!work||count<n*n||!finite_array(a,n*n))return TINYLINALG_INVALID_ARGUMENT;
    memcpy(work,a,n*n*sizeof *work);
    for(size_t k=0;k<n;++k){
        size_t pivot=k;double best=fabs(work[k*n+k]);
        for(size_t i=k+1;i<n;++i)if(fabs(work[i*n+k])>best){best=fabs(work[i*n+k]);pivot=i;}
        if(best==0)return TINYLINALG_SINGULAR;
        p[k]=pivot;
        if(pivot!=k)for(size_t j=0;j<n;++j){double x=work[k*n+j];work[k*n+j]=work[pivot*n+j];work[pivot*n+j]=x;}
        for(size_t i=k+1;i<n;++i){
            work[i*n+k]/=work[k*n+k];
            for(size_t j=k+1;j<n;++j)work[i*n+j]-=work[i*n+k]*work[k*n+j];
        }
        if(!finite_array(work,n*n))return TINYLINALG_ARITHMETIC_ERROR;
    }
    memmove(lu,work,n*n*sizeof *lu);memcpy(pivots,p,n*sizeof *pivots);return TINYLINALG_OK;
}
tinylinalg_status tinylinalg_lu_solve(size_t n,const double *lu,const size_t *pivots,
    const double *b,double *x,double *work,size_t count) {
    if(!valid_size(n)||!pivots||!x||!work||count<n||!finite_array(lu,n*n)||!finite_array(b,n))return TINYLINALG_INVALID_ARGUMENT;
    for(size_t k=0;k<n;++k)if(pivots[k]<k||pivots[k]>=n)return TINYLINALG_INVALID_ARGUMENT;
    memcpy(work,b,n*sizeof *work);
    for(size_t k=0;k<n;++k){double z=work[k];work[k]=work[pivots[k]];work[pivots[k]]=z;}
    for(size_t i=0;i<n;++i)for(size_t j=0;j<i;++j)work[i]-=lu[i*n+j]*work[j];
    for(size_t ii=n;ii>0;--ii){
        size_t i=ii-1;
        if(lu[i*n+i]==0)return TINYLINALG_SINGULAR;
        for(size_t j=i+1;j<n;++j)work[i]-=lu[i*n+j]*work[j];
        work[i]/=lu[i*n+i];
    }
    if(!finite_array(work,n))return TINYLINALG_ARITHMETIC_ERROR;
    memcpy(x,work,n*sizeof *x);return TINYLINALG_OK;
}
tinylinalg_status tinylinalg_solve(size_t n,const double *a,const double *b,
    double *x,double *work,size_t count) {
    size_t pivots[TINYLINALG_MAX_DIM];
    tinylinalg_status status;
    if(!valid_size(n)||!work||count<n*n+n||!x||!finite_array(b,n))return TINYLINALG_INVALID_ARGUMENT;
    status=tinylinalg_lu(n,a,work,pivots,work,count);
    if(status!=TINYLINALG_OK)return status;
    return tinylinalg_lu_solve(n,work,pivots,b,x,work+n*n,count-n*n);
}
static int symmetric(size_t n,const double *a,double *scale) {
    *scale=0;
    for(size_t i=0;i<n*n;++i)*scale=fmax(*scale,fabs(a[i]));
    for(size_t i=0;i<n;++i)for(size_t j=i+1;j<n;++j)
        if(fabs(a[i*n+j]-a[j*n+i])>64*DBL_EPSILON*(*scale))return 0;
    return 1;
}
tinylinalg_status tinylinalg_cholesky(size_t n,const double *a,double *lower,
    double *work,size_t count) {
    double scale;
    if(!valid_size(n)||!lower||!work||count<n*n||!finite_array(a,n*n)||!symmetric(n,a,&scale))return TINYLINALG_INVALID_ARGUMENT;
    for(size_t i=0;i<n*n;++i)work[i]=0;
    for(size_t i=0;i<n;++i)for(size_t j=0;j<=i;++j){
        double value=a[i*n+j];
        for(size_t k=0;k<j;++k)value-=work[i*n+k]*work[j*n+k];
        if(j==i){if(value<=0)return TINYLINALG_NOT_POSITIVE_DEFINITE;work[i*n+j]=sqrt(value);}
        else work[i*n+j]=value/work[j*n+j];
        if(!isfinite(work[i*n+j]))return TINYLINALG_ARITHMETIC_ERROR;
    }
    memcpy(lower,work,n*n*sizeof *lower);return TINYLINALG_OK;
}
tinylinalg_status tinylinalg_eigen_symmetric(size_t n,const double *a,
    double *values,double *vectors,double *work,size_t count) {
    double scale,evalues[TINYLINALG_MAX_DIM];
    double *m=work,*v;
    int converged=0;
    if(!valid_size(n)||!values||!vectors||!work||count<2*n*n||!finite_array(a,n*n)||!symmetric(n,a,&scale))return TINYLINALG_INVALID_ARGUMENT;
    v=work+n*n;
    /* Scale before iteration to avoid overflow in 2*off-diagonal and differences. */
    if(scale==0)scale=1;
    for(size_t i=0;i<n;++i)for(size_t j=0;j<n;++j){m[i*n+j]=a[i*n+j]/scale;v[i*n+j]=(i==j)?1:0;}
    for(size_t sweep=0;sweep<64;++sweep){
        double off=0;
        for(size_t p=0;p<n;++p)for(size_t q=p+1;q<n;++q)off=fmax(off,fabs(m[p*n+q]));
        if(off<=16*DBL_EPSILON){converged=1;break;}
        for(size_t p=0;p<n;++p)for(size_t q=p+1;q<n;++q){
            double apq=m[p*n+q],tau,t,c,s,app,aqq;
            if(fabs(apq)<=DBL_EPSILON)continue;
            tau=(m[q*n+q]-m[p*n+p])/(2*apq);
            t=copysign(1.0,tau)/(fabs(tau)+hypot(1.0,tau));c=1/sqrt(1+t*t);s=t*c;
            app=m[p*n+p];aqq=m[q*n+q];
            m[p*n+p]=app-t*apq;m[q*n+q]=aqq+t*apq;m[p*n+q]=m[q*n+p]=0;
            for(size_t k=0;k<n;++k)if(k!=p&&k!=q){
                double mkp=m[k*n+p],mkq=m[k*n+q];
                m[k*n+p]=m[p*n+k]=c*mkp-s*mkq;m[k*n+q]=m[q*n+k]=s*mkp+c*mkq;
            }
            for(size_t k=0;k<n;++k){double vp=v[k*n+p],vq=v[k*n+q];v[k*n+p]=c*vp-s*vq;v[k*n+q]=s*vp+c*vq;}
        }
    }
    if(!converged)return TINYLINALG_NO_CONVERGENCE;
    for(size_t i=0;i<n;++i){evalues[i]=m[i*n+i]*scale;if(!isfinite(evalues[i]))return TINYLINALG_ARITHMETIC_ERROR;}
    for(size_t i=0;i<n;++i){
        size_t k=i;for(size_t j=i+1;j<n;++j)if(evalues[j]<evalues[k])k=j;
        if(k!=i){double x=evalues[i];evalues[i]=evalues[k];evalues[k]=x;for(size_t row=0;row<n;++row){x=v[row*n+i];v[row*n+i]=v[row*n+k];v[row*n+k]=x;}}
    }
    memcpy(values,evalues,n*sizeof *values);memcpy(vectors,v,n*n*sizeof *vectors);return TINYLINALG_OK;
}
