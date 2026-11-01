/* SPDX-License-Identifier: MIT */
#include "tinysgp4fit/tinysgp4fit.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "data/reference.h"
#define CHECK(x) do { if (!(x)) {fprintf(stderr,"line %d: %s\n",__LINE__,#x);return 1;} } while(0)
int main(void)
{
    tinysgp4fit_workspace workspace;
    double max_position=0,max_velocity=0,max_covariance=0,max_cost=0;
    size_t i,j,k;
    for(i=0;i<sizeof reference/sizeof reference[0];++i) {
        tinysgp4fit_options options;tinysgp4fit_result result;tinysgp4fit_status status;
        double actual[7],expected[7],tolerances[7]={2e-7,2e-7,2e-9,2e-7,2e-7,1e-9,2e-8};
        CHECK(tinysgp4fit_default_options(&options)==TINYSGP4FIT_OK);
        options.fit_bstar=reference[i].fit_bstar;options.use_velocity=reference[i].use_velocity;
        status=tinysgp4fit_fit(&reference[i].seed,reference[i].observations,8,&options,&workspace,&result);
        if(status!=TINYSGP4FIT_OK)fprintf(stderr,"fit case %zu status %d\n",i,(int)status);
        CHECK(status==TINYSGP4FIT_OK);
        actual[0]=result.elements.inclination_deg;actual[1]=result.elements.raan_deg;actual[2]=result.elements.eccentricity;actual[3]=result.elements.arg_perigee_deg;actual[4]=result.elements.mean_anomaly_deg;actual[5]=result.elements.mean_motion_rev_day;actual[6]=result.elements.bstar;
        expected[0]=reference[i].fitted.inclination_deg;expected[1]=reference[i].fitted.raan_deg;expected[2]=reference[i].fitted.eccentricity;expected[3]=reference[i].fitted.arg_perigee_deg;expected[4]=reference[i].fitted.mean_anomaly_deg;expected[5]=reference[i].fitted.mean_motion_rev_day;expected[6]=reference[i].fitted.bstar;
        for(j=0;j<7;++j){if(fabs(actual[j]-expected[j])>=tolerances[j])fprintf(stderr,"case%zu param%zu actual%.17g expected%.17g error%.9g\n",i,j,actual[j],expected[j],fabs(actual[j]-expected[j]));CHECK(fabs(actual[j]-expected[j])<tolerances[j]);}
        CHECK(result.elements.catalog_id==reference[i].seed.catalog_id);
        CHECK(result.elements.epoch_jd1==reference[i].seed.epoch_jd1&&result.elements.epoch_jd2==reference[i].seed.epoch_jd2);
        CHECK(result.elements.ephemeris_type==reference[i].seed.ephemeris_type);
        CHECK(result.parameter_count==(options.fit_bstar?7u:6u));
        CHECK(result.scalar_observation_count==(options.use_velocity?48u:24u));
        if(!options.fit_bstar)CHECK(result.elements.bstar==reference[i].seed.bstar);
        CHECK(result.weighted_sum_squares<reference[i].seed_cost*1e-6);
        CHECK(fabs(result.elements.inclination_deg-reference[i].seed.inclination_deg)>.001);
        {double e=fabs(result.weighted_sum_squares-reference[i].cost);if(e>max_cost)max_cost=e;CHECK(e<1e-5);}
        CHECK(fabs(result.normalized_rms-sqrt(result.weighted_sum_squares/result.scalar_observation_count))<1e-14);
        for(j=0;j<result.parameter_count;++j)for(k=0;k<result.parameter_count;++k){size_t n=result.parameter_count;double scale=sqrt(fabs(reference[i].covariance[j*n+j]))*sqrt(fabs(reference[i].covariance[k*n+k]));double error=fabs(result.covariance[j*n+k]-reference[i].covariance[j*n+k])/fmax(scale,1e-30);if(error>max_covariance)max_covariance=error;if(error>=1e-4)fprintf(stderr,"cov case%zu j%zu k%zu %.9g\n",i,j,k,error);CHECK(error<1e-4);}
        CHECK(tinytle_initialize_elements(&result.elements,&workspace.satellite)==TINYTLE_OK);
        for(j=0;j<4;++j){CHECK(tgp4_propagate(&workspace.satellite,reference[i].holdout[j][0],&workspace.state)==TGP4_OK);for(k=0;k<3;++k){double er=fabs(workspace.state.position_km[k]*1000-reference[i].holdout[j][1+k]);double ev=fabs(workspace.state.velocity_km_s[k]*1000-reference[i].holdout[j][4+k]);if(er>max_position)max_position=er;if(ev>max_velocity)max_velocity=ev;if(er>=.001||ev>=1e-6)fprintf(stderr,"holdout case%zu %.9g m %.9g m/s\n",i,er,ev);CHECK(er<.001);CHECK(ev<1e-6);}}
    }
    {
        tinysgp4fit_options options;tinysgp4fit_result out={0},saved;tinysgp4fit_observation bad[8];tinytle_elements seed=reference[0].seed;
        memcpy(&saved,&out,sizeof saved);CHECK(tinysgp4fit_default_options(&options)==TINYSGP4FIT_OK);
        memcpy(bad,observations,sizeof bad);bad[0].position_sigma_m[0]=0;
        CHECK(tinysgp4fit_fit(&seed,bad,8,&options,&workspace,&out)==TINYSGP4FIT_INVALID);CHECK(!memcmp(&out,&saved,sizeof out));
        memcpy(bad,observations,sizeof bad);for(i=0;i<8;++i)bad[i].minutes_since_epoch=0;options.fit_bstar=1;
        CHECK(tinysgp4fit_fit(&seed,bad,8,&options,&workspace,&out)==TINYSGP4FIT_UNOBSERVABLE);CHECK(!memcmp(&out,&saved,sizeof out));
        options.fit_bstar=0;seed.eccentricity=0;CHECK(tinysgp4fit_fit(&seed,observations,8,&options,&workspace,&out)==TINYSGP4FIT_UNOBSERVABLE);CHECK(!memcmp(&out,&saved,sizeof out));
        seed=reference[0].seed;options.max_iterations=1;
        CHECK(tinysgp4fit_fit(&seed,observations,8,&options,&workspace,&out)==TINYSGP4FIT_NO_CONVERGENCE);CHECK(!memcmp(&out,&saved,sizeof out));
    }
    printf("SGP4 fit:7 SciPy/python-sgp4 fits,28 holdouts; max %.9g m / %.9g m/s; normalized covariance %.9g; cost %.9g; workspace %zu bytes\n",max_position,max_velocity,max_covariance,max_cost,sizeof workspace);
    return 0;
}
