/* SPDX-License-Identifier: MIT */
#include "tinypkepler/tinypkepler.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "data/reference.h"
#define CHECK(x) do { if (!(x)) {fprintf(stderr,"line %d: %s\n",__LINE__,#x);return 1;} } while(0)
int main(void)
{
    size_t i,j;double max_r=0,max_v=0;
    for(i=0;i<sizeof reference/sizeof reference[0];++i){const double *r=reference[i];tinyelements_classical e={r[3],r[4],r[5],r[6],r[7],r[8]};tinypkepler_result out;CHECK(tinypkepler_propagate(r[0],r[1],r[2],&e,r[9],&out)==TINYPKEPLER_OK);for(j=0;j<3;++j){double er=fabs(out.state.position_m[j]-r[10+j]),ev=fabs(out.state.velocity_m_s[j]-r[13+j]);if(er>max_r)max_r=er;if(ev>max_v)max_v=ev;if(er>=.001||ev>=1e-6)fprintf(stderr,"case %zu %.9g m %.9g m/s\n",i,er,ev);CHECK(er<.001);CHECK(ev<1e-6);}CHECK(out.mean_elements.p_m==e.p_m && out.mean_elements.eccentricity==e.eccentricity && out.mean_elements.inclination_rad==e.inclination_rad);}
    {tinyelements_classical e={7000000,1,1,0,0,0};tinypkepler_result out={0},saved=out;CHECK(tinypkepler_propagate(3.986e14,6378137,.001,&e,60,&out)==TINYPKEPLER_INVALID);CHECK(!memcmp(&out,&saved,sizeof out));e.eccentricity=.1;CHECK(tinypkepler_propagate(3.986e14,6378137,0,&e,60,&out)==TINYPKEPLER_OK);CHECK(out.raan_rate_rad_s==0&&out.arg_periapsis_rate_rad_s==0);}
    printf("Perturbed Kepler: 36 pinned Vallado J2 cases; max %.9g m / %.9g m/s\n",max_r,max_v);return 0;
}
