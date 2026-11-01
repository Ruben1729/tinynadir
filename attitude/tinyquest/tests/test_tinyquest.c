/* SPDX-License-Identifier: MIT */
#include "tinyquest/tinyquest.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL %s:%d: %s\n",__FILE__,__LINE__,#c);return 1;}}while(0)
#include "data/reference.h"
static double qdistance(const double a[4],const double b[4]){double m=0,p=0;for(size_t i=0;i<4;++i){m+=(a[i]-b[i])*(a[i]-b[i]);p+=(a[i]+b[i])*(a[i]+b[i]);}return sqrt(fmin(m,p));}
int main(void){
    tinyquest_workspace work;double q[4],loss,max_q=0,max_loss=0;
    for(size_t c=0;c<sizeof reference_cases/sizeof reference_cases[0];++c){
        CHECK(tinyquest_solve(reference_cases[c].count,reference_cases[c].reference,reference_cases[c].body,reference_cases[c].weights,q,&loss,&work)==TINYQUEST_OK);
        double error=qdistance(q,reference_cases[c].quaternion);CHECK(error<3e-12);max_q=fmax(max_q,error);
        CHECK(fabs(loss-reference_cases[c].loss)<3e-13);max_loss=fmax(max_loss,fabs(loss-reference_cases[c].loss));
    }
    {double r[]={1,0,0,0,1,0},b[]={0,1,0,-1,0,0},w[]={DBL_MAX,DBL_MAX};
      CHECK(tinyquest_solve(2,r,b,w,q,&loss,&work)==TINYQUEST_OK);CHECK(fabs(q[0]-sqrt(.5))<1e-13&&fabs(q[3]-sqrt(.5))<1e-13);
      CHECK(tinyquest_solve(2,r,b,NULL,q,NULL,&work)==TINYQUEST_OK);}
    {double r[]={1,0,0,1,0,0},b[]={0,1,0,0,1,0},w[]={1,1},old[4]={123,456,789,42};memcpy(q,old,sizeof q);loss=123;
      CHECK(tinyquest_solve(2,r,b,w,q,&loss,&work)==TINYQUEST_DEGENERATE);CHECK(!memcmp(q,old,sizeof q)&&loss==123);
      w[0]=w[1]=0;CHECK(tinyquest_solve(2,r,b,w,q,&loss,&work)==TINYQUEST_DEGENERATE);
      w[0]=-1;CHECK(tinyquest_solve(2,r,b,w,q,&loss,&work)==TINYQUEST_INVALID_ARGUMENT);
      w[0]=NAN;CHECK(tinyquest_solve(2,r,b,w,q,&loss,&work)==TINYQUEST_INVALID_ARGUMENT);
      CHECK(tinyquest_solve(1,r,b,NULL,q,&loss,&work)==TINYQUEST_INVALID_ARGUMENT);
      CHECK(tinyquest_solve(65,r,b,NULL,q,&loss,&work)==TINYQUEST_INVALID_ARGUMENT);
      CHECK(tinyquest_solve(2,NULL,b,NULL,q,&loss,&work)==TINYQUEST_INVALID_ARGUMENT);
      CHECK(tinyquest_solve(2,r,b,NULL,q,&loss,NULL)==TINYQUEST_INVALID_ARGUMENT);}
    printf("tinyquest: 96 SciPy Kabsch cases; max quaternion difference %.3g, normalized loss %.3g\n",max_q,max_loss);return 0;
}
