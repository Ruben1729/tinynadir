#ifndef TGP4_H
#define TGP4_H
#include <stddef.h>
#include "internal.h"
#include "status.h"
/* Internal record layout is not a stable API. No shared mutable state or heap. */
typedef struct { tgp4_record record; int initialized; } tgp4_satellite;
typedef struct { double position_km[3]; double velocity_km_s[3]; } tgp4_state;
/* Zero-initialize satellite before first use. 69-byte TLE lines, optional CR/LF.
   Numeric IDs, ephemeris type 0, WGS-72, improved operation mode.
   Failed calls preserve the satellite and output. Different satellites can be
   used concurrently; serialize calls that access the same satellite. */
#ifdef __cplusplus
extern "C" {
#endif
tgp4_status tgp4_initialize(tgp4_satellite*,const char*,size_t,const char*,size_t);
/* Minutes relative to the TLE epoch. Finite requests within +/-10080 only.
   Outputs are TEME position (km) and velocity (km/s). */
tgp4_status tgp4_propagate(tgp4_satellite*,double minutes_since_epoch,tgp4_state*);
#ifdef __cplusplus
}
#endif
#endif
