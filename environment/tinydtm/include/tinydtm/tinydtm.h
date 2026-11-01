#ifndef TINYDTM_H
#define TINYDTM_H
#ifdef __cplusplus
extern "C" {
#endif
/* Selected model is DTM2000, NOT DTM2012/DTM2020. SI/radian interface. */
typedef struct {
 int day_of_year;
 double altitude_m,longitude_rad,latitude_rad,local_solar_time_rad;
 double f107,f107_mean,kp_3h,kp_24h;
} tinydtm_input;
typedef struct {
 double density_kg_m3,temperature_k,exospheric_temperature_k;
 /* Order H,He,O,N2,O2,N. */
 double number_density_m3[6],partial_mass_density_kg_m3[6];
} tinydtm_output;
/* Per-call workspace, layout is implementation detail; no initialization required. */
typedef struct {
 int day;double f[3],fbar[3],akp[5];
 double clfl,slfl,p10,p20,p30,p40,p50,p60,p11,p21,p31,p41,p51;
 double p22,p32,p42,p52,p62,p33,p10mg,p20mg,p40mg,hl0,ch,sh,c2h,s2h,c3h,s3h;
 double scratch[97];
} tinydtm_workspace;
typedef enum {TINYDTM_OK=0,TINYDTM_INVALID,TINYDTM_RANGE,TINYDTM_ARITHMETIC} tinydtm_status;
/* Profile altitude120..1500 km, F10.7 40..400, Kp0..9. No hidden data lookup. */
tinydtm_status tinydtm_evaluate(const tinydtm_input *input,tinydtm_workspace *workspace,tinydtm_output *out);
#ifdef __cplusplus
}
#endif
#endif
