/* SPDX-License-Identifier: MIT; adapted original model has public-domain notice. */
#include "tinymsis/tinymsis.h"
#include "core.h"
#include <math.h>
#include <string.h>
#include <stddef.h>
void tinymsis_default_switches(int switches[24]) {int i;if(switches)for(i=0;i<24;i++)switches[i]=1;}
tinymsis_status tinymsis_evaluate(const tinymsis_input *in,const int switches[24],tinymsis_workspace *ctx,tinymsis_output *out) {
    struct nrlmsise_input model;struct nrlmsise_flags flags;struct nrlmsise_output result;struct ap_array ap;
    tinymsis_output r;int i;static const int index[8]={0,1,2,3,4,6,7,8};
    if(!in||!ctx||!out)return TINYMSIS_INVALID;
    if(!isfinite(in->seconds_utc)||!isfinite(in->altitude_m)||!isfinite(in->latitude_rad)||!isfinite(in->longitude_rad)||!isfinite(in->local_solar_time_hours)||!isfinite(in->f107_mean)||!isfinite(in->f107))return TINYMSIS_INVALID;
    for(i=0;i<7;i++)if(!isfinite(in->ap[i]))return TINYMSIS_INVALID;
    if(in->day_of_year<1||in->day_of_year>366||in->seconds_utc<0||in->seconds_utc>=86400||in->altitude_m<0||in->altitude_m>1500000||fabs(in->latitude_rad)>1.57079632679489661923||in->local_solar_time_hours<0||in->local_solar_time_hours>24||in->f107_mean<40||in->f107_mean>400||in->f107<40||in->f107>400||(in->effective_drag_density!=0&&in->effective_drag_density!=1))return TINYMSIS_RANGE;
    for(i=0;i<7;i++)if(in->ap[i]<0||in->ap[i]>400)return TINYMSIS_RANGE;
    for(i=0;i<24;i++) {
        int s=switches?switches[i]:1;
        if((i==0&&s!=1)||(i!=0&&(s<0||s>2)&&!(i==9&&s==-1)))return TINYMSIS_RANGE;
        flags.switches[i]=s;
    }
    memset(ctx,0,sizeof *ctx);memset(&result,0,sizeof result);
    model.year=0;model.doy=in->day_of_year;model.sec=in->seconds_utc;model.alt=in->altitude_m/1000;
    model.g_lat=in->latitude_rad*(180/3.141592653589793238462643383279502884);
    model.g_long=remainder(in->longitude_rad,6.28318530717958647693)*(180/3.141592653589793238462643383279502884);
    model.lst=in->local_solar_time_hours;model.f107A=in->f107_mean;model.f107=in->f107;model.ap=in->ap[0];
    for(i=0;i<7;i++)ap.a[i]=in->ap[i];model.ap_a=&ap;
    if(in->effective_drag_density)tinymsis_core_gtd7d(ctx,&model,&flags,&result);else tinymsis_core_gtd7(ctx,&model,&flags,&result);
    if(ctx->numeric_error)return TINYMSIS_ARITHMETIC;
    for(i=0;i<8;i++){r.number_density_m3[i]=result.d[index[i]];if(!isfinite(r.number_density_m3[i])||r.number_density_m3[i]<0)return TINYMSIS_ARITHMETIC;}
    r.density_kg_m3=result.d[5];r.exospheric_temperature_k=result.t[0];r.temperature_k=result.t[1];
    if(!isfinite(r.density_kg_m3)||r.density_kg_m3<0||!isfinite(r.exospheric_temperature_k)||r.exospheric_temperature_k<=0||!isfinite(r.temperature_k)||r.temperature_k<=0)return TINYMSIS_ARITHMETIC;
    *out=r;return TINYMSIS_OK;
}
