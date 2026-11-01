#include "tinymeanosculating/tinymeanosculating.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define PI 3.14159265358979323846
#define CHECK(x) do{if(!(x)){fprintf(stderr,"%s:%d %s\n",__FILE__,__LINE__,#x);exit(1);}}while(0)
#include "data/reference.h"
static void pack(const tinymeanosculating_mean*q,double x[6]){x[0]=q->semimajor_axis_m;x[1]=q->f;x[2]=q->g;x[3]=q->h;x[4]=q->k;x[5]=q->mean_longitude_rad;}
static void unpack(const double*x,tinymeanosculating_mean*q){q->semimajor_axis_m=x[0];q->f=x[1];q->g=x[2];q->h=x[3];q->k=x[4];q->mean_longitude_rad=x[5];}
int main(void){
 tinymeanosculating_workspace w;tinymeanosculating_options o={128,2.,1000000,12,2e-12};tinymeanosculating_model m;tinycowell_state s,back,savedstate;tinymeanosculating_mean avg,reference,saved;tinymeanosculating_stats stats,savedstats;size_t i,maxsteps=0;unsigned maxpasses=0;int k;double maxa=0,maxangle=0,maxposition=0,maxvelocity=0;
 for(i=0;i<sizeof mean_cases/sizeof mean_cases[0];i++){
  const double*q=mean_cases[i];double actual[6];m.mu_m3_s2=q[0];m.equatorial_radius_m=q[1];m.j2=q[2];o.samples=(unsigned)q[9];
  for(k=0;k<3;k++){s.position_m[k]=q[3+k];s.velocity_m_s[k]=q[6+k];}
  {tinymeanosculating_status code=tinymeanosculating_from_osculating(&m,&s,&o,&w,&avg,&stats);if(code)fprintf(stderr,"case %zu status %d\n",i,(int)code);CHECK(code==TINYMEANOSCULATING_OK);}
  pack(&avg,actual);
  for(k=0;k<6;k++){
   double diff=fabs(k==5?remainder(actual[k]-q[10+k],2*PI):actual[k]-q[10+k]);
   if(k==0){CHECK(diff<.001);if(diff>maxa)maxa=diff;}else {CHECK(diff<1e-10);if(diff>maxangle)maxangle=diff;}
  }
  CHECK(stats.averaging_passes==1&&stats.integration_steps<=o.max_steps);
  unpack(q+10,&reference);
  CHECK(tinymeanosculating_to_osculating(&m,&reference,&o,&w,&back,&stats)==TINYMEANOSCULATING_OK);
  if(stats.integration_steps>maxsteps)maxsteps=stats.integration_steps;if(stats.averaging_passes>maxpasses)maxpasses=stats.averaging_passes;
  CHECK(stats.residual<=o.inverse_tolerance);
  for(k=0;k<3;k++){
   double pe=fabs(back.position_m[k]-s.position_m[k]),ve=fabs(back.velocity_m_s[k]-s.velocity_m_s[k]);
   CHECK(pe<.002);CHECK(ve<2e-6);if(pe>maxposition)maxposition=pe;if(ve>maxvelocity)maxvelocity=ve;
  }
  CHECK(tinymeanosculating_from_osculating(&m,&back,&o,&w,&avg,&stats)==TINYMEANOSCULATING_OK);pack(&avg,actual);
  CHECK(fabs(actual[0]-reference.semimajor_axis_m)/reference.semimajor_axis_m<=o.inverse_tolerance);
 }
 /* All error paths must preserve result/stats, including exhausted budgets/inverse. */
 saved=avg;savedstats=stats;o.samples=127;
 CHECK(tinymeanosculating_from_osculating(&m,&s,&o,&w,&avg,&stats)==TINYMEANOSCULATING_RANGE);CHECK(!memcmp(&avg,&saved,sizeof avg)&&!memcmp(&stats,&savedstats,sizeof stats));
 o.samples=128;o.max_steps=1;CHECK(tinymeanosculating_from_osculating(&m,&s,&o,&w,&avg,&stats)==TINYMEANOSCULATING_BUDGET);CHECK(!memcmp(&avg,&saved,sizeof avg)&&!memcmp(&stats,&savedstats,sizeof stats));
 o.max_steps=1000000;m.j2=.00108262668;o.max_iterations=1;savedstate=back;
 CHECK(tinymeanosculating_to_osculating(&m,&reference,&o,&w,&back,&stats)==TINYMEANOSCULATING_NO_CONVERGENCE);CHECK(!memcmp(&back,&savedstate,sizeof back)&&!memcmp(&stats,&savedstats,sizeof stats));
 o.max_iterations=12;s.position_m[0]=NAN;CHECK(tinymeanosculating_from_osculating(&m,&s,&o,&w,&avg,&stats)==TINYMEANOSCULATING_INVALID);
 printf("tinymeanosculating: %zu independent SciPy averages + inversions; max %.9g m a, %.9g coordinate, %.9g m position, %.9g m/s; max %zu steps/%u passes\n",sizeof mean_cases/sizeof mean_cases[0],maxa,maxangle,maxposition,maxvelocity,maxsteps,maxpasses);return 0;
}
