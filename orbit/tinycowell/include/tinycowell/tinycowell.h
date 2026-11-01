#ifndef TINYCOWELL_H
#define TINYCOWELL_H
#include <stddef.h>
#include "tinyrkf45/tinyrkf45.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct {double position_m[3],velocity_m_s[3];} tinycowell_state;
/* Caller-supplied additional acceleration in the same inertial frame. */
typedef int (*tinycowell_force)(double epoch_s,const tinycowell_state *state,void *context,double acceleration_m_s2[3]);
typedef struct {double mu_m3_s2;tinycowell_force force;void *context;} tinycowell_model;
typedef enum {TINYCOWELL_OK=0,TINYCOWELL_INVALID,TINYCOWELL_FORCE_ERROR,TINYCOWELL_INTEGRATION_ERROR} tinycowell_status;
/* mu>=0 adds central gravity; use mu=0 if callback already supplies the full
   acceleration (e.g. tinygravity). Coordinates must be inertial. No Coriolis
   or rotating-frame terms are implied. Outputs preserved on failure.
   No events/maneuvers. Workspace must be disjoint from all state/output. */
tinycowell_status tinycowell_fixed(const tinycowell_model *model,double start_s,double end_s,const tinycowell_state *initial,double max_step_s,size_t max_steps,tinycowell_state *out,double *workspace,size_t count);
/* Fixed path workspace48 doubles, adaptive54 doubles; initial/out may alias. */
tinycowell_status tinycowell_adaptive(const tinycowell_model *model,double start_s,double end_s,const tinycowell_state *initial,const tinyrkf45_options *options,tinycowell_state *out,tinyrkf45_stats *stats,double *workspace,size_t count);
#ifdef __cplusplus
}
#endif
#endif
