/* SPDX-License-Identifier: MIT */
#include "tinytle/tinytle.h"
tinytle_status tinytle_initialize(const char *a,size_t na,const char *b,size_t nb,tgp4_satellite *sat){tinytle_elements g;tinytle_status status=tinytle_parse(a,na,b,nb,&g);if(status!=TINYTLE_OK)return status;return tinytle_initialize_elements(&g,sat);}
