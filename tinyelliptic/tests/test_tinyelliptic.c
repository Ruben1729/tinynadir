#include "tinyelliptic/tinyelliptic.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL %d: %s\n",__LINE__,#c);return 1;}}while(0)
#include "data/reference.h"
int main(void){double value,maximum=0;
 for(size_t i=0;i<sizeof integral_cases/sizeof integral_cases[0];++i){tinyelliptic_status s=tinyelliptic_first(integral_cases[i].phi,integral_cases[i].m,&value);if(integral_cases[i].singular)CHECK(s==TINYELLIPTIC_SINGULAR);else{CHECK(s==TINYELLIPTIC_OK);double error=fabs(value-integral_cases[i].first);CHECK(error<8e-14*fmax(1,fabs(integral_cases[i].first)));maximum=fmax(maximum,error);}CHECK(tinyelliptic_second(integral_cases[i].phi,integral_cases[i].m,&value)==TINYELLIPTIC_OK);CHECK(fabs(value-integral_cases[i].second)<8e-14*fmax(1,fabs(integral_cases[i].second)));}
 for(size_t i=0;i<sizeof complete_cases/sizeof complete_cases[0];++i){if(complete_cases[i].m==1)CHECK(tinyelliptic_complete_first(1,&value)==TINYELLIPTIC_SINGULAR);else{CHECK(tinyelliptic_complete_first(complete_cases[i].m,&value)==TINYELLIPTIC_OK);CHECK(fabs(value-complete_cases[i].first)<8e-14*fmax(1,complete_cases[i].first));}CHECK(tinyelliptic_complete_second(complete_cases[i].m,&value)==TINYELLIPTIC_OK);CHECK(fabs(value-complete_cases[i].second)<8e-14);}
 for(size_t i=0;i<sizeof inverse_cases/sizeof inverse_cases[0];++i){CHECK(tinyelliptic_inverse_second(inverse_cases[i].value,inverse_cases[i].m,&value)==TINYELLIPTIC_OK);CHECK(fabs(value-inverse_cases[i].phi)<5e-12);double recovered;CHECK(tinyelliptic_second(value,inverse_cases[i].m,&recovered)==TINYELLIPTIC_OK);CHECK(fabs(recovered-inverse_cases[i].value)<5e-13);}
 for(size_t i=0;i<sizeof arc_cases/sizeof arc_cases[0];++i){CHECK(tinyelliptic_arc(arc_cases[i].a,arc_cases[i].b,arc_cases[i].start,arc_cases[i].end,&value)==TINYELLIPTIC_OK);CHECK(fabs(value-arc_cases[i].length)<5e-12*fmax(1,fabs(arc_cases[i].length)));}
 {value=123;CHECK(tinyelliptic_first(1,1.1,&value)==TINYELLIPTIC_INVALID_ARGUMENT&&value==123);CHECK(tinyelliptic_second(NAN,.5,&value)==TINYELLIPTIC_INVALID_ARGUMENT);CHECK(tinyelliptic_inverse_second(1,-.1,&value)==TINYELLIPTIC_INVALID_ARGUMENT);CHECK(tinyelliptic_arc(0,1,0,1,&value)==TINYELLIPTIC_INVALID_ARGUMENT);CHECK(tinyelliptic_first(1e7,.1,&value)==TINYELLIPTIC_RANGE);CHECK(tinyelliptic_arc(2,2,0,3,&value)==TINYELLIPTIC_OK&&value==6);}
 printf("tinyelliptic: 210 incomplete/11 complete/49 inverse/64 arc SciPy cases; max first-integral difference %.3g\n",maximum);return 0;}
