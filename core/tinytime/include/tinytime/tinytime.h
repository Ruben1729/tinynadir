#ifndef TINYTIME_H
#define TINYTIME_H
#ifdef __cplusplus
extern "C" {
#endif
/* Uniform 86400-second days. first is midnight Julian date (... .5),
   fraction is in [0,1). The time scale is selected by the caller/API. */
typedef struct { double first, fraction; } tinytime_jd;
typedef enum { TINYTIME_OK=0, TINYTIME_INVALID, TINYTIME_RANGE } tinytime_status;
tinytime_status tinytime_normalize(double first, double second, tinytime_jd *out);
tinytime_status tinytime_calendar_to_jd(int year,int month,int day,int hour,int minute,double second,tinytime_jd *out);
tinytime_status tinytime_jd_to_calendar(tinytime_jd jd,int *year,int *month,int *day,int *hour,int *minute,double *second);
tinytime_status tinytime_add_seconds(tinytime_jd jd,double seconds,tinytime_jd *out);
tinytime_status tinytime_difference_seconds(tinytime_jd a,tinytime_jd b,double *seconds);
tinytime_status tinytime_decimal_year(tinytime_jd jd,double *year);
/* UTC coordinate conversions require caller-supplied TAI-UTC and UT1-UTC.
   This uniform-day API rejects second=60; use a leap-aware clock adapter at
   an actual UTC leap-second instant. No internal leap-second table. */
tinytime_status tinytime_utc_to_tai(tinytime_jd utc,double tai_minus_utc_s,tinytime_jd *tai);
tinytime_status tinytime_utc_to_ut1(tinytime_jd utc,double ut1_minus_utc_s,tinytime_jd *ut1);
tinytime_status tinytime_tai_to_tt(tinytime_jd tai,tinytime_jd *tt);
tinytime_status tinytime_tt_to_tai(tinytime_jd tt,tinytime_jd *tai);
/* Geocentric two-term periodic approximation, not a relativistic/topocentric
   clock model. Validated 1900-2100; |TDB-TT error| < 0.1 ms in sampled cases. */
tinytime_status tinytime_tt_to_tdb(tinytime_jd tt,tinytime_jd *tdb);
#ifdef __cplusplus
}
#endif
#endif
