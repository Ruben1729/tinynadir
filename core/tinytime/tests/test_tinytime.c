#include "tinytime/tinytime.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(x) do {if(!(x)){fprintf(stderr,"%s:%d: %s\n",__FILE__,__LINE__,#x);exit(1);}}while(0)
#include "data/calendar.h"
#include "data/tdb.h"
int main(void) {
    size_t i;double max_tdb=0,delta; tinytime_jd q,z,saved;
    for(i=0;i<sizeof calendar_cases/sizeof calendar_cases[0];i++) {
        int y,m,d,h,min;double s,decimal;
        CHECK(tinytime_calendar_to_jd(calendar_cases[i].y,calendar_cases[i].m,calendar_cases[i].d,calendar_cases[i].h,calendar_cases[i].min,calendar_cases[i].s,&q)==TINYTIME_OK);
        CHECK(q.first==calendar_cases[i].first);CHECK(fabs(q.fraction-calendar_cases[i].fraction)<2e-16);
        CHECK(tinytime_jd_to_calendar(q,&y,&m,&d,&h,&min,&s)==TINYTIME_OK);
        CHECK(y==calendar_cases[i].y && m==calendar_cases[i].m && d==calendar_cases[i].d && h==calendar_cases[i].h && min==calendar_cases[i].min);
        CHECK(fabs(s-calendar_cases[i].s)<2e-11);
        CHECK(tinytime_decimal_year(q,&decimal)==TINYTIME_OK);CHECK(decimal>=y && decimal<=y+1.0);
        CHECK(tinytime_add_seconds(q,12345.678,&z)==TINYTIME_OK);CHECK(tinytime_difference_seconds(z,q,&delta)==TINYTIME_OK);CHECK(fabs(delta-12345.678)<2e-10);
    }
    CHECK(tinytime_calendar_to_jd(2000,1,1,12,0,0,&q)==TINYTIME_OK);
    CHECK(tinytime_utc_to_tai(q,37,&z)==TINYTIME_OK);CHECK(tinytime_difference_seconds(z,q,&delta)==TINYTIME_OK);CHECK(fabs(delta-37)<1e-10);
    CHECK(tinytime_tai_to_tt(z,&z)==TINYTIME_OK);CHECK(tinytime_difference_seconds(z,q,&delta)==TINYTIME_OK);CHECK(fabs(delta-69.184)<1e-10);
    CHECK(tinytime_tt_to_tai(z,&z)==TINYTIME_OK);CHECK(tinytime_difference_seconds(z,q,&delta)==TINYTIME_OK);CHECK(fabs(delta-37)<1e-10);
    CHECK(tinytime_utc_to_ut1(q,-.234,&z)==TINYTIME_OK);CHECK(tinytime_difference_seconds(z,q,&delta)==TINYTIME_OK);CHECK(fabs(delta+.234)<1e-10);
    for(i=0;i<sizeof tdb_cases/sizeof tdb_cases[0];i++) {
        double e; q.first=tdb_cases[i][0];q.fraction=tdb_cases[i][1];
        CHECK(tinytime_tt_to_tdb(q,&z)==TINYTIME_OK);CHECK(tinytime_difference_seconds(z,q,&delta)==TINYTIME_OK);
        e=fabs(delta-tdb_cases[i][2]);CHECK(e<0.0001);if(e>max_tdb)max_tdb=e;
    }
    q.first=123;q.fraction=.123;saved=q;
    CHECK(tinytime_calendar_to_jd(1900,2,29,0,0,0,&q)==TINYTIME_RANGE);CHECK(!memcmp(&q,&saved,sizeof q));
    CHECK(tinytime_calendar_to_jd(2025,1,1,0,0,60,&q)==TINYTIME_RANGE);
    CHECK(tinytime_calendar_to_jd(2000,1,1,0,0,NAN,&q)==TINYTIME_INVALID);
    CHECK(tinytime_normalize(INFINITY,0,&q)==TINYTIME_INVALID);CHECK(!memcmp(&q,&saved,sizeof q));
    CHECK(tinytime_add_seconds(q,NAN,&q)==TINYTIME_INVALID);
    CHECK(tinytime_calendar_to_jd(2026,11,1,0,0,0,NULL)==TINYTIME_INVALID);
    printf("tinytime: ERFA calendar/scales and frozen TDB comparisons pass; max sampled TDB error %.9g seconds\n",max_tdb);
    return 0;
}
