#ifndef TINYOMM_H
#define TINYOMM_H
#include <stddef.h>
#include "tinytle/tinytle.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef enum {TINYOMM_OK=0,TINYOMM_INVALID,TINYOMM_UNSUPPORTED,TINYOMM_PROPAGATION_ERROR} tinyomm_status;
typedef struct {const char *name;size_t name_length;const char *value;size_t value_length;} tinyomm_field;
/* Decoded-field adapter or bounded KVN subset. Keys are uppercase CCSDS/GP
   names. Fields may reference nonterminated caller storage. Max64 fields,
   each name<=64 and value<=128; KVN<=8192bytes/line<=256bytes.
   Required: EPOCH, MEAN_MOTION, ECCENTRICITY, INCLINATION, RA_OF_ASC_NODE,
   ARG_OF_PERICENTER, MEAN_ANOMALY, BSTAR. Strict mode also requires
   CENTER_NAME=EARTH, REF_FRAME=TEME, TIME_SYSTEM=UTC, MEAN_ELEMENT_THEORY=SGP4.
   celestrak_defaults=1 explicitly supplies those four GP endpoint defaults.
   Other supplied incompatible metadata/type/version is always rejected.
   Unknown metadata ignored, duplicate known keys rejected. Type0 only.
   ISO epochs have4digit Gregorian years and optional fractional seconds/Z;
   actual leap seconds(second60) are not accepted by the uniform-day profile.
   Outputs preserved on failure; raw JSON/XML/CSV decoding is outside this API. */
tinyomm_status tinyomm_parse_fields(const tinyomm_field *fields,size_t count,int celestrak_defaults,tinytle_elements *out);
tinyomm_status tinyomm_parse_kvn(const char *text,size_t length,int celestrak_defaults,tinytle_elements *out);
tinyomm_status tinyomm_initialize_kvn(const char *text,size_t length,int celestrak_defaults,tgp4_satellite *satellite);
#ifdef __cplusplus
}
#endif
#endif
