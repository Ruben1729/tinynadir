/* SPDX-License-Identifier: MIT */
#include "tinysgp4/tinysgp4.h"
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL %s:%d: %s\n",__FILE__,__LINE__,#c); exit(1); } } while (0)
typedef struct {
    char a[70],b[70];
    tgp4_status initialization;
    /* B*, eccentricity, inclination, node, argument, mean anomaly, mean motion. */
    double elements[7];
} fixture;
#include "data/cases.h"
#define FIXTURE_COUNT (sizeof fixtures/sizeof fixtures[0])
#define epoch_issue fixtures[33]

static void checksum(char* p) {
    int sum=0;
    for(int i=0;i<68;++i) { if(p[i]>='0'&&p[i]<='9')sum+=p[i]-'0';else if(p[i]=='-')++sum; }
    p[68]=(char)('0'+sum%10);p[69]=0;
}
static tgp4_status init(tgp4_satellite* satellite,const fixture* f) {
    return tgp4_initialize(satellite,f->a,69,f->b,69);
}
static double distance(const double a[3],const double b[3]) {
    double sum=0;for(int i=0;i<3;++i)sum+=(a[i]-b[i])*(a[i]-b[i]);return sqrt(sum);
}
static void replay(const char* path,char mode,unsigned expected_states,unsigned expected_errors) {
    /* These are test-only records, not storage required by the library. Multiple
       records preserve the original interleaved request/cache regression. */
    tgp4_satellite satellites[FIXTURE_COUNT];
    memset(satellites,0,sizeof satellites);
    unsigned init_errors=0;
    for(size_t i=0;i<FIXTURE_COUNT;++i) {
        CHECK(init(&satellites[i],&fixtures[i])==fixtures[i].initialization);
        if(fixtures[i].initialization!=TGP4_OK) {
            tgp4_satellite empty;memset(&empty,0,sizeof empty);CHECK(!memcmp(&satellites[i],&empty,sizeof empty));
            ++init_errors;continue;
        }
        const tgp4_record* r=&satellites[i].record;
        const double actual[]={r->bstar,r->ecco,r->inclo,r->nodeo,r->argpo,r->mo,r->no_kozai};
        for(size_t k=0;k<7;++k)CHECK(actual[k]==fixtures[i].elements[k]);
        /* The older published vectors use legacy mode. Public API defaults to
           improved mode; only this explicit compatibility test selects legacy. */
        if(mode=='a')CHECK(tgp4_core_sgp4init('a',&satellites[i].record));
    }
    CHECK(init_errors==1);
    FILE* input=fopen(path,"r");CHECK(input);
    unsigned states=0,errors=0;double maximum_position=0,maximum_velocity=0;
    char line[512];
    while(fgets(line,sizeof line,input)) {
        unsigned index;int reference_error;double t,r[3],v[3];
        CHECK(strchr(line,'\n')!=NULL);
        CHECK(sscanf(line,"%u,%lf,%d,%lf,%lf,%lf,%lf,%lf,%lf",&index,&t,&reference_error,
                     &r[0],&r[1],&r[2],&v[0],&v[1],&v[2])==9);
        CHECK(index<FIXTURE_COUNT);CHECK(fixtures[index].initialization==TGP4_OK);
        CHECK(isfinite(t));CHECK(t>=-10080 && t<=10080);
        tgp4_satellite before;memcpy(&before,&satellites[index],sizeof before);
        tgp4_state state={{123,456,789},{42,43,44}},unchanged=state;
        tgp4_status status=tgp4_propagate(&satellites[index],t,&state);
        if(reference_error) {
            CHECK(status==TGP4_PROPAGATION_ERROR);
            CHECK(!memcmp(&state,&unchanged,sizeof state));
            CHECK(!memcmp(&satellites[index],&before,sizeof before));++errors;continue;
        }
        CHECK(status==TGP4_OK);
        double er=distance(state.position_km,r)*1000,ev=distance(state.velocity_km_s,v)*1000;
        CHECK(isfinite(er)&&isfinite(ev));
        /* Extreme synthetic coordinates can make one binary64 ULP exceed a
           millimetre. Keep the absolute floor plus a small relative allowance. */
        double position_limit=0.001+32*DBL_EPSILON*sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])*1000;
        double velocity_limit=0.00001+32*DBL_EPSILON*sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])*1000;
        if(er>=position_limit || ev>=velocity_limit)
            fprintf(stderr,"fixture %u t=%.17g: %.12g m / %.12g m/s\n",index,t,er,ev);
        CHECK(er<position_limit);CHECK(ev<velocity_limit);
        if(er>maximum_position)maximum_position=er;
        if(ev>maximum_velocity)maximum_velocity=ev;++states;
    }
    CHECK(!ferror(input));fclose(input);
    CHECK(states==expected_states);CHECK(errors==expected_errors);
    printf("%s: %u matching states, %u expected errors; max %.9g m / %.9g m/s\n",
           path,states,errors,maximum_position,maximum_velocity);
}
static void epoch_regression(void) {
    /* vaelen/sgp4 PR #1, opened 2026-09-01: formatter adds half a day twice.
       This library has no formatter. Assert the exact report's parsed epoch. */
    tgp4_satellite sat={0};CHECK(init(&sat,&epoch_issue)==TGP4_OK);
    CHECK(sat.record.jdsatepoch==2461283.5);
    CHECK(fabs(sat.record.jdsatepochF-0.16371918)<1e-13);
    double noon_minutes=((2461283.5-sat.record.jdsatepoch)+(0.5-sat.record.jdsatepochF))*1440;
    CHECK(fabs(noon_minutes-484.2443808)<1e-9);
    /* Independent Gregorian anchors cover TLE century pivot, leap years,
       midnight/noon, end of day, and end-of-year normalization. */
    const struct {const char* epoch;double jd,fraction;int valid;} dates[]={
        {"00001.00000000",2451544.5,0,1},
        {"00060.50000000",2451603.5,0.5,1},
        {"00366.99999999",2451909.5,0.99999999,1},
        {"01366.00000000",0,0,0},
        {"00000.50000000",0,0,0},
        {"00367.00000000",0,0,0},
        {"57001.00000000",2435839.5,0,1},
        {"56001.00000000",2471998.5,0,1},
        {"26243.49999999",2461283.5,0.49999999,1},
        {"26243.50000000",2461283.5,0.5,1},
        {"26243.99999999",2461283.5,0.99999999,1}
    };
    for(size_t i=0;i<sizeof dates/sizeof dates[0];++i) {
        fixture f=epoch_issue;memcpy(f.a+18,dates[i].epoch,14);checksum(f.a);
        tgp4_status status=init(&sat,&f);
        if(!dates[i].valid){CHECK(status==TGP4_INVALID_TLE);continue;}
        CHECK(status==TGP4_OK);CHECK(sat.record.jdsatepoch==dates[i].jd);
        CHECK(fabs(sat.record.jdsatepochF-dates[i].fraction)<1e-13);
    }
    puts("Epoch regression: exact September 2026 report and Gregorian boundary cases pass");
}
static void parser_failures(void) {
    tgp4_satellite sat={0};CHECK(init(&sat,&epoch_issue)==TGP4_OK);
    tgp4_satellite saved;memcpy(&saved,&sat,sizeof saved);
    for(size_t n=0;n<69;++n) {
        CHECK(tgp4_initialize(&sat,epoch_issue.a,n,epoch_issue.b,69)==TGP4_INVALID_TLE);
        CHECK(tgp4_initialize(&sat,epoch_issue.a,69,epoch_issue.b,n)==TGP4_INVALID_TLE);
    }
    CHECK(tgp4_initialize(NULL,epoch_issue.a,69,epoch_issue.b,69)==TGP4_INVALID_TLE);
    CHECK(tgp4_initialize(&sat,NULL,69,epoch_issue.b,69)==TGP4_INVALID_TLE);
    CHECK(tgp4_initialize(&sat,epoch_issue.a,69,NULL,69)==TGP4_INVALID_TLE);
    fixture bad=epoch_issue;bad.b[2]='9';checksum(bad.b);CHECK(init(&sat,&bad)==TGP4_INVALID_TLE);
    bad=epoch_issue;memcpy(bad.a+18,"-1",2);checksum(bad.a);CHECK(init(&sat,&bad)==TGP4_INVALID_TLE);
    bad=epoch_issue;memcpy(bad.b+8,"181.0000",8);checksum(bad.b);CHECK(init(&sat,&bad)==TGP4_INVALID_TLE);
    bad=epoch_issue;memcpy(bad.b+52," 0.00000000",11);checksum(bad.b);CHECK(init(&sat,&bad)==TGP4_INVALID_TLE);
    bad=epoch_issue;bad.a[62]='4';checksum(bad.a);CHECK(init(&sat,&bad)==TGP4_UNSUPPORTED_EPHEMERIS);
    for(int i=0;i<69;++i) {
        bad=epoch_issue;bad.a[i]=(bad.a[i]=='9'?'8':'9');
        CHECK(init(&sat,&bad)==TGP4_INVALID_TLE);
    }
    CHECK(!memcmp(&sat,&saved,sizeof sat));
    char raw_a[69],raw_b[69];memcpy(raw_a,epoch_issue.a,69);memcpy(raw_b,epoch_issue.b,69);
    CHECK(tgp4_initialize(&sat,raw_a,sizeof raw_a,raw_b,sizeof raw_b)==TGP4_OK);
    char a[71],b[71];memcpy(a,epoch_issue.a,69);memcpy(b,epoch_issue.b,69);
    a[69]=b[69]='\r';a[70]=b[70]='\n';
    CHECK(tgp4_initialize(&sat,a,71,b,71)==TGP4_OK);
    a[70]='X';CHECK(tgp4_initialize(&sat,a,71,b,71)==TGP4_INVALID_TLE);
    tgp4_state out={{1,2,3},{4,5,6}},before=out;
    CHECK(tgp4_propagate(NULL,0,&out)==TGP4_UNINITIALIZED);
    CHECK(tgp4_propagate(&sat,0,NULL)==TGP4_UNINITIALIZED);
    CHECK(tgp4_propagate(&sat,NAN,&out)==TGP4_TIME_OUT_OF_RANGE);
    CHECK(tgp4_propagate(&sat,-INFINITY,&out)==TGP4_TIME_OUT_OF_RANGE);
    CHECK(tgp4_propagate(&sat,-10080.1,&out)==TGP4_TIME_OUT_OF_RANGE);
    CHECK(!memcmp(&out,&before,sizeof out));
    puts("Parser/API: short and nonterminated buffers, checksums, nulls, fields, CRLF, failure preservation pass");
}
int main(void) {
    CHECK(FIXTURE_COUNT==34);
    replay("tests/data/vallado.csv",'i',19530,1400);
    replay("tests/data/published.csv",'a',597,0);
    epoch_regression();parser_failures();
    printf("Storage on this ABI: satellite=%zu bytes, output=%zu bytes\n",sizeof(tgp4_satellite),sizeof(tgp4_state));
    return 0;
}
