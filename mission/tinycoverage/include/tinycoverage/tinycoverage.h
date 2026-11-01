#ifndef TINYCOVERAGE_H
#define TINYCOVERAGE_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYCOVERAGE_OK=0,TINYCOVERAGE_INVALID,TINYCOVERAGE_RANGE,TINYCOVERAGE_GEOMETRY,TINYCOVERAGE_CAPACITY,TINYCOVERAGE_BUDGET,TINYCOVERAGE_CALLBACK,TINYCOVERAGE_RESOLUTION} tinycoverage_status;
typedef struct {double surface_distance_m,forward_azimuth_rad;} tinycoverage_arc;
typedef struct {double latitude_rad,longitude_rad;} tinycoverage_point;
typedef struct {double central_angle_rad,surface_radius_m,area_m2;int horizon_limited;} tinycoverage_footprint;
typedef struct {double range_m,azimuth_rad,elevation_rad;} tinycoverage_observation;
typedef struct {double nodal_period_s,semimajor_axis_m;} tinycoverage_repeat;
tinycoverage_status tinycoverage_great_circle(double radius_m,double latitude1_rad,double longitude1_rad,double latitude2_rad,double longitude2_rad,tinycoverage_arc *out);
tinycoverage_status tinycoverage_destination(double radius_m,double latitude_rad,double longitude_rad,double azimuth_rad,double surface_distance_m,tinycoverage_point *out);
tinycoverage_status tinycoverage_nadir_footprint(double radius_m,double altitude_m,double cone_half_angle_rad,tinycoverage_footprint *out);
/* Station position is ECEF; local geodetic basis is supplied independently. */
tinycoverage_status tinycoverage_observe(const double station_ecef_m[3],double station_latitude_rad,double station_longitude_rad,const double object_ecef_m[3],tinycoverage_observation *out);
/* Ideal nodal repeat: revolutions periods = rotations effective Earth days.
   Use an effective rotation period accounting for node drift if appropriate.
   Semi-major axis estimate is two-body and excludes J2. */
tinycoverage_status tinycoverage_repeat_orbit(double mu_m3_s2,double effective_rotation_period_s,unsigned revolutions,unsigned rotations,tinycoverage_repeat *out);
#define TINYCOVERAGE_MAX_INTERVALS 64
typedef struct {double entry_s,exit_s;} tinycoverage_interval;
typedef struct {tinycoverage_interval intervals[TINYCOVERAGE_MAX_INTERVALS];} tinycoverage_workspace;
/* Callback returns0 on success and writes finite scalar margin; margin>=0 is
   inside access (e.g. elevation minus minimum elevation). No hidden propagator. */
typedef int (*tinycoverage_margin_fn)(double seconds,void *context,double *margin);
/* Fixed-step sign-change search plus bisection. Events shorter than sample_s
   can be missed; it is the application's responsibility to choose sampling.
   count/calls/output intervals are committed only on complete success.
   scratch may change on any call. Bisection capped at64 iterations/crossing. */
tinycoverage_status tinycoverage_scan(tinycoverage_margin_fn margin,void *context,double start_s,double end_s,double sample_s,double time_tolerance_s,size_t max_calls,tinycoverage_workspace *scratch,tinycoverage_interval *intervals,size_t capacity,size_t *count,size_t *calls);
#ifdef __cplusplus
}
#endif
#endif
