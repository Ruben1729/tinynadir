/* SPDX-License-Identifier: MIT */
#include "tinycovariance/tinycovariance.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "data/reference.h"
#define CHECK(x) do { if (!(x)) {fprintf(stderr,"line %d: %s\n",__LINE__,#x);return 1;} } while(0)
static double scaled_error(const double *a,const double *b,size_t n)
{size_t i,j;double max=0;for(i=0;i<n;++i)for(j=0;j<n;++j){double scale=sqrt(fabs(b[i*n+i]))*sqrt(fabs(b[j*n+j]));double e=fabs(a[i*n+j]-b[i*n+j])/fmax(scale,1e-30);if(e>max)max=e;}return max;}
static int nonlinear(const double *x,size_t n,double *y,size_t m,void *context)
{(void)context;if(n!=2||m!=2)return 0;y[0]=x[0]*x[0]+x[1];y[1]=atan2(x[1],x[0]);return 1;}
int main(void)
{
    size_t i,j;double maximum=0,out[36];
    CHECK(tinycovariance_transform(3,2,generic_J,generic_P,out)==TINYCOVARIANCE_OK);CHECK(scaled_error(out,generic_output,2)<1e-14);
    CHECK(tinycovariance_rotate_state(rotation,cartesian_P,out)==TINYCOVARIANCE_OK);CHECK(scaled_error(out,rotated_P,6)<1e-13);
    for(i=0;i<sizeof reference/sizeof reference[0];++i){tinyelements_state state,back;double P[36],again[36],error;for(j=0;j<3;++j){state.position_m[j]=reference[i].state[j];state.velocity_m_s[j]=reference[i].state[3+j];}if(reference[i].mode==0){tinyelements_classical e;CHECK(tinycovariance_cartesian_to_classical(reference[i].mu,&state,reference[i].P,&e,P)==TINYCOVARIANCE_OK);CHECK(fabs(e.p_m-reference[i].elements[0])<1e-6);CHECK(tinycovariance_classical_to_cartesian(reference[i].mu,&e,P,&back,again)==TINYCOVARIANCE_OK);}else{tinyelements_equinoctial e;CHECK(tinycovariance_cartesian_to_equinoctial(reference[i].mu,&state,reference[i].P,&e,P)==TINYCOVARIANCE_OK);CHECK(fabs(e.p_m-reference[i].elements[0])<1e-6);CHECK(tinycovariance_equinoctial_to_cartesian(reference[i].mu,&e,P,&back,again)==TINYCOVARIANCE_OK);}error=scaled_error(P,reference[i].covariance,6);if(error>maximum)maximum=error;if(error>=3e-6)fprintf(stderr,"case %zu scaled cov error %.9g\n",i,error);CHECK(error<3e-6);CHECK(scaled_error(again,reference[i].P,6)<1e-5);}
    {double x[2]={2,1},step[2]={1e-4,1e-4},period[2]={0,6.2831853071795864769},J[4];CHECK(tinycovariance_jacobian(2,2,x,step,period,nonlinear,NULL,J)==TINYCOVARIANCE_OK);CHECK(fabs(J[0]-4)<1e-10);CHECK(fabs(J[1]-1)<1e-10);CHECK(fabs(J[2]+.2)<1e-9);CHECK(fabs(J[3]-.4)<1e-9);}
    {double P[4]={1,2,2,1},J[4]={1,0,0,1},save[4]={9,8,7,6};for(i=0;i<4;++i)out[i]=save[i];CHECK(tinycovariance_transform(2,2,J,P,out)==TINYCOVARIANCE_NOT_PSD);CHECK(!memcmp(save,out,sizeof save));}
    printf("Covariance: NumPy congruence/rotation and6 SciPy/Vallado orbital maps; max normalized error %.9g\n",maximum);return 0;
}
