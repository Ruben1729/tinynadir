#ifndef TINYTLE_H
#define TINYTLE_H
#include <stddef.h>
#include <stdint.h>
#include "tinysgp4/tinysgp4.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYTLE_OK=0,TINYTLE_INVALID,TINYTLE_UNSUPPORTED,TINYTLE_PROPAGATION_ERROR} tinytle_status;
/* Canonical conventional-GP input fields. Format-native angular units are
   degrees and revolutions/day, not SI orbit-state units. Epoch is UTC with
   midnight Julian date ending .5 and fraction in[0,1). B* is inverse Earth
   radii; catalog id is metadata and not used in orbital equations. */
typedef struct {uint32_t catalog_id;int ephemeris_type;double epoch_jd1,epoch_jd2,inclination_deg,raan_deg,eccentricity,arg_perigee_deg,mean_anomaly_deg,mean_motion_rev_day,bstar;} tinytle_elements;
/* Exactly69 columns plus optional CR/LF. Buffers need not be terminated.
   Supports numeric0..99999 (leading zeros/spaces) and Alpha-5 A0000..Z9999
   excluding I/O. Blank BSTAR is zero; epoch day000 is the prior Dec31.
   Checksums, matching ids, field ranges, ephemeris type0. Failed out preserved. */
tinytle_status tinytle_parse(const char *line1,size_t length1,const char *line2,size_t length2,tinytle_elements *out);
/* SGP4 initialization from decoded fields, WGS72 improved mode, preserving
   Vallado unit-conversion evaluation order. This is a separate archive object
   from the text parser. Decoded epochs1600..9999. Successful initialization
   validates a finite epoch TEME state. Failed satellite preserved; no allocations. */
tinytle_status tinytle_initialize_elements(const tinytle_elements *elements,tgp4_satellite *satellite);
tinytle_status tinytle_initialize(const char *line1,size_t length1,const char *line2,size_t length2,tgp4_satellite *satellite);
#ifdef __cplusplus
}
#endif
#endif
