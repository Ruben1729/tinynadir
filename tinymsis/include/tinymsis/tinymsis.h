#ifndef TINYMSIS_H
#define TINYMSIS_H
#ifdef __cplusplus
extern "C" {
#endif
/* Selected atmosphere model is NRLMSISE-00, not restricted NRLMSIS 2.x. */
typedef struct {
 int day_of_year;
 double seconds_utc,altitude_m,latitude_rad,longitude_rad,local_solar_time_hours;
 double f107_mean,f107;
 /* daily Ap, current3h, previous3h/6h/9h, means12..33h and36..57h */
 double ap[7];
 int effective_drag_density;
} tinymsis_input;
typedef struct {
 /* He,O,N2,O2,Ar,H,N,anomalous O; number per cubic metre. */
 double number_density_m3[8];
 double density_kg_m3,exospheric_temperature_k,temperature_k;
} tinymsis_output;
/* Caller scratch; layout is an implementation detail. Reset each evaluation.
   Separate workspaces permit independent/reentrant satellite calls. */
typedef struct {
    double gsurf;
    double re;
    double dd;
    double dm04, dm16, dm28, dm32, dm40, dm01, dm14;
    double meso_tn1[5];
    double meso_tn2[4];
    double meso_tn3[5];
    double meso_tgn1[2];
    double meso_tgn2[2];
    double meso_tgn3[2];
    double dfa;
    double plg[4][9];
    double ctloc, stloc;
    double c2tloc, s2tloc;
    double s3tloc, c3tloc;
    double apdf, apt[4];
    int numeric_error;
} tinymsis_workspace;
typedef enum {TINYMSIS_OK=0,TINYMSIS_INVALID,TINYMSIS_RANGE,TINYMSIS_ARITHMETIC} tinymsis_status;
/* Default switches[0]=1 (SI), other switches1. Optional switch9=-1 uses Ap history. */
void tinymsis_default_switches(int switches[24]);
/* switches may be NULL for defaults. Otherwise switch0 must remain1; switches
   1..23 accept0/1/2, with -1 additionally allowed for switch9. Altitude0..1500km. */
tinymsis_status tinymsis_evaluate(const tinymsis_input *input,const int switches[24],tinymsis_workspace *workspace,tinymsis_output *out);
#ifdef __cplusplus
}
#endif
#endif
