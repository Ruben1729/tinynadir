/* SPDX-License-Identifier: MIT */
#include "tinyomm/tinyomm.h"
tinyomm_status tinyomm_initialize_kvn(const char *text,size_t n,int defaults,tgp4_satellite *sat){tinytle_elements g;tinytle_status init;tinyomm_status status=tinyomm_parse_kvn(text,n,defaults,&g);if(status!=TINYOMM_OK)return status;init=tinytle_initialize_elements(&g,sat);return init==TINYTLE_OK?TINYOMM_OK:(init==TINYTLE_UNSUPPORTED?TINYOMM_UNSUPPORTED:(init==TINYTLE_INVALID?TINYOMM_INVALID:TINYOMM_PROPAGATION_ERROR));}
