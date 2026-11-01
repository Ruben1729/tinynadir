/* Copyright 2002-2026 CS GROUP
 * Licensed to CS GROUP (CS) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * CS licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/* TinyNadir C99 adaptation of scalar DTM2000 computation only; dynamic arrays
   replaced by caller-owned bounded workspace and immutable coefficients. */
#include "tinydtm/tinydtm.h"
#include <math.h>
#include <stddef.h>
#include "coefficients.h"
#define RE 6356.77
#define ZLB0 120.0
#define CPMG .19081
#define SPMG .98163
#define XLMG -1.2392
#define GSURF 980.665
#define RGAS 831.4
#define ROT .017214206
#define ROT2 .034428412
static const double ALEFA[]={0,-.40,-.38,0,0,0,0};
static const double MA[]={0,1,4,16,28,32,14};
static const double VMA[]={0,1.6606e-24,6.6423e-24,26.569e-24,46.4958e-24,53.1381e-24,23.2479e-24};

static double semestrial(tinydtm_workspace *ctx,double param) {

            int debeq_pr = 59;
            int debeq_au = 244;
            double result;
            if (ctx->day >= 100) {
                double xmult  = (ctx->day - debeq_au) / 40.0;
                result = param - 2.0 * param * xmult;
            } else {
                double xmult  = (ctx->day - debeq_pr) / 40.0;
                result = 2.0 * param * xmult - param;
            }
            return result;

}

static double g_function(tinydtm_workspace *ctx,const double a[97],double da[97],int ff0,int kle_eq) {


            double fmfb[3]={0};
            double fbm150[3]={0};

            // latitude terms
            da[2]  = ctx->p20;
            da[3]  = ctx->p40;
            da[74] = ctx->p10;
            double a74 = a[74];
            double a77 = a[77];
            double a78 = a[78];
            if (kle_eq == -1) {
                // winter
                a74 = -a74;
                a77 = -a77;
                a78 = -a78;
            }
            if (kle_eq == 0 ) {
                // equinox
                a74 = semestrial(ctx,a74);
                a77 = semestrial(ctx,a77);
                a78 = semestrial(ctx,a78);
            }
            da[77] = ctx->p30;
            da[78] = ctx->p50;
            da[79] = ctx->p60;

            // flux terms
            fmfb[1]   = ctx->f[1] - ctx->fbar[1];
            fmfb[2]   = ctx->f[2] - ctx->fbar[2];
            fbm150[1] = ctx->fbar[1] - 150.0;
            fbm150[2] = ctx->fbar[2];
            da[4]     = fmfb[1];
            da[6]     = fbm150[1];
            da[4]     = da[4] + a[70] * fmfb[2];
            da[6]     = da[6] + a[71] * fbm150[2];
            da[70]    = fmfb[2] * (a[4] + 2.0 * a[5] * da[4] + a[82] * ctx->p10 +
                                   a[83] * ctx->p20 + a[84] * ctx->p30);
            da[71]    = fbm150[2] * (a[6] + 2.0 * a[69] * da[6] + a[85] * ctx->p10 +
                                     a[86] * ctx->p20 + a[87] * ctx->p30);
            da[5]     = da[4] * da[4];
            da[69]    = da[6] * da[6];
            da[82]    = da[4] * ctx->p10;
            da[83]    = da[4] * ctx->p20;
            da[84]    = da[4] * ctx->p30;
            da[85]    = da[6] * ctx->p20;
            da[86]    = da[6] * ctx->p30;
            da[87]    = da[6] * ctx->p40;

            // Kp terms
            int ikp  = 62;
            int ikpm = 67;
            double c2fi = 1.0 - ctx->p10mg * ctx->p10mg;
            double dkp  = ctx->akp[1] + (a[ikp] + c2fi * a[ikp + 1]) * ctx->akp[2];
            double dakp = a[7] + a[8] * ctx->p20mg + a[68] * ctx->p40mg +
                          2.0 * dkp * (a[60] + a[61] * ctx->p20mg +
                                       a[75] * 2.0 * dkp * dkp);
            da[ikp] = dakp * ctx->akp[2];
            da[ikp + 1] = da[ikp] * c2fi;
            double dkpm  = ctx->akp[3] + a[ikpm] * ctx->akp[4];
            double dakpm = a[64] + a[65] * ctx->p20mg + a[72] * ctx->p40mg +
                                 2.0 * dkpm * (a[66] + a[73] * ctx->p20mg +
                                               a[76] * 2.0 * dkpm * dkpm);
            da[ikpm] = dakpm * ctx->akp[4];
            da[7]    = dkp;
            da[8]    = ctx->p20mg * dkp;
            da[68]   = ctx->p40mg * dkp;
            da[60]   = dkp * dkp;
            da[61]   = ctx->p20mg * da[60];
            da[75]   = da[60] * da[60];
            da[64]   = dkpm;
            da[65]   = ctx->p20mg * dkpm;
            da[72]   = ctx->p40mg * dkpm;
            da[66]   = dkpm * dkpm;
            da[73]   = ctx->p20mg * da[66];
            da[76]   = da[66] * da[66];

            // non-periodic g(l) function
            double f0 = a[4]  * da[4]  + a[5]  * da[5]  + a[6]  * da[6]  +
                        a[69] * da[69] + a[82] * da[82] + a[83] * da[83] +
                        a[84] * da[84] + a[85] * da[85] + a[86] * da[86] +
                        a[87] * da[87];
            double f1f = 1.0 + f0 * ff0;

            f0 = f0 + a[2] * da[2] + a[3] * da[3] + a74 * da[74] +
                 a77 * da[77] + a[7] * da[7] + a[8] * da[8] +
                 a[60] * da[60] + a[61] * da[61] + a[68] * da[68] +
                 a[64] * da[64] + a[65] * da[65] + a[66] * da[66] +
                 a[72] * da[72] + a[73] * da[73] + a[75] * da[75] +
                 a[76] * da[76] + a78   * da[78] + a[79] * da[79];
//          termes annuels symetriques en latitude
            da[9]  = cos(ROT * (ctx->day - a[11]));
            da[10] = ctx->p20 * da[9];
//          termes semi-annuels symetriques en latitude
            da[12] = cos(ROT2 * (ctx->day - a[14]));
            da[13] = ctx->p20 * da[12];
//          termes annuels non symetriques en latitude
            double coste = cos(ROT * (ctx->day - a[18]));
            da[15] = ctx->p10 * coste;
            da[16] = ctx->p30 * coste;
            da[17] = ctx->p50 * coste;
//          terme  semi-annuel  non symetrique  en latitude
            double cos2te = cos(ROT2 * (ctx->day - a[20]));
            da[19] = ctx->p10 * cos2te;
            da[39] = ctx->p30 * cos2te;
            da[59] = ctx->p50 * cos2te;
//          termes diurnes [et couples annuel]
            da[21] = ctx->p11 * ctx->ch;
            da[22] = ctx->p31 * ctx->ch;
            da[23] = ctx->p51 * ctx->ch;
            da[24] = da[21] * coste;
            da[25] = ctx->p21 * ctx->ch * coste;
            da[26] = ctx->p11 * ctx->sh;
            da[27] = ctx->p31 * ctx->sh;
            da[28] = ctx->p51 * ctx->sh;
            da[29] = da[26] * coste;
            da[30] = ctx->p21 * ctx->sh * coste;
//          termes semi-diurnes [et couples annuel]
            da[31] = ctx->p22 * ctx->c2h;
            da[37] = ctx->p42 * ctx->c2h;
            da[32] = ctx->p32 * ctx->c2h * coste;
            da[33] = ctx->p22 * ctx->s2h;
            da[38] = ctx->p42 * ctx->s2h;
            da[34] = ctx->p32 * ctx->s2h * coste;
            da[88] = ctx->p32 * ctx->c2h;
            da[89] = ctx->p32 * ctx->s2h;
            da[90] = ctx->p52 * ctx->c2h;
            da[91] = ctx->p52 * ctx->s2h;
            double a88 = a[88];
            double a89 = a[89];
            double a90 = a[90];
            double a91 = a[91];
            if (kle_eq == -1) {            //hiver
                a88 = -a88;
                a89 = -a89;
                a90 = -a90;
                a91 = -a91;
            }
            if (kle_eq == 0) {             //equinox
                a88 = semestrial(ctx,a88);
                a89 = semestrial(ctx,a89);
                a90 = semestrial(ctx,a90);
                a91 = semestrial(ctx,a91);
            }
            da[92] = ctx->p62 * ctx->c2h;
            da[93] = ctx->p62 * ctx->s2h;
//          termes ter-diurnes
            da[35] = ctx->p33 * ctx->c3h;
            da[36] = ctx->p33 * ctx->s3h;
//          fonction g[l] periodique
            double fp = a[9]  * da[9]  + a[10] * da[10] + a[12] * da[12] + a[13] * da[13] +
                        a[15] * da[15] + a[16] * da[16] + a[17] * da[17] + a[19] * da[19] +
                        a[21] * da[21] + a[22] * da[22] + a[23] * da[23] + a[24] * da[24] +
                        a[25] * da[25] + a[26] * da[26] + a[27] * da[27] + a[28] * da[28] +
                        a[29] * da[29] + a[30] * da[30] + a[31] * da[31] + a[32] * da[32] +
                        a[33] * da[33] + a[34] * da[34] + a[35] * da[35] + a[36] * da[36] +
                        a[37] * da[37] + a[38] * da[38] + a[39] * da[39] + a[59] * da[59] +
                        a88   * da[88] + a89   * da[89] + a90   * da[90] + a91   * da[91] +
                        a[92] * da[92] + a[93] * da[93];
//          termes d'activite magnetique
            da[40] = ctx->p10 * coste * dkp;
            da[41] = ctx->p30 * coste * dkp;
            da[42] = ctx->p50 * coste * dkp;
            da[43] = ctx->p11 * ctx->ch * dkp;
            da[44] = ctx->p31 * ctx->ch * dkp;
            da[45] = ctx->p51 * ctx->ch * dkp;
            da[46] = ctx->p11 * ctx->sh * dkp;
            da[47] = ctx->p31 * ctx->sh * dkp;
            da[48] = ctx->p51 * ctx->sh * dkp;

//          fonction g[l] periodique supplementaire
            fp += a[40] * da[40] + a[41] * da[41] + a[42] * da[42] + a[43] * da[43] +
                  a[44] * da[44] + a[45] * da[45] + a[46] * da[46] + a[47] * da[47] +
                  a[48] * da[48];

            dakp = (a[40] * ctx->p10 + a[41] * ctx->p30 + a[42] * ctx->p50) * coste +
                   (a[43] * ctx->p11 + a[44] * ctx->p31 + a[45] * ctx->p51) * ctx->ch +
                   (a[46] * ctx->p11 + a[47] * ctx->p31 + a[48] * ctx->p51) * ctx->sh;
            da[ikp] += dakp * ctx->akp[2];
            da[ikp + 1] = da[ikp] + dakp * c2fi * ctx->akp[2];
//          termes de longitude
            da[49] = ctx->p11 * ctx->clfl;
            da[50] = ctx->p21 * ctx->clfl;
            da[51] = ctx->p31 * ctx->clfl;
            da[52] = ctx->p41 * ctx->clfl;
            da[53] = ctx->p51 * ctx->clfl;
            da[54] = ctx->p11 * ctx->slfl;
            da[55] = ctx->p21 * ctx->slfl;
            da[56] = ctx->p31 * ctx->slfl;
            da[57] = ctx->p41 * ctx->slfl;
            da[58] = ctx->p51 * ctx->slfl;

//          fonction g[l] periodique supplementaire
            fp += a[49] * da[49] + a[50] * da[50] + a[51] * da[51] + a[52] * da[52] +
                  a[53] * da[53] + a[54] * da[54] + a[55] * da[55] + a[56] * da[56] +
                  a[57] * da[57] + a[58] * da[58];

//          fonction g(l) totale (couplage avec le flux)
            return f0 + fp * f1f;


}

static void compute(tinydtm_workspace *ctx,double altiKM,double lon,double lat,double hl,tinydtm_output *out) {
    double *da=ctx->scratch;







            // Sine and cosine of local latitude and longitude
            double scLat_s=sin(lat),scLat_c=cos(lat);
            double scLon_s=sin(lon),scLon_c=cos(lon);

            // compute Legendre polynomials wrt geographic pole
            double c = scLat_s;
            double c2 = c * c;
            double c4 = c2 * c2;
            double s = scLat_c;
            double s2 = s * s;
            ctx->p10 = c;
            ctx->p20 = 1.5 * c2 - 0.5;
            ctx->p30 = c * (2.5 * c2 - 1.5);
            ctx->p40 = 4.375 * c4 - 3.75 * c2 + 0.375;
            ctx->p50 = c * (7.875 * c4 - 8.75 * c2 + 1.875);
            ctx->p60 = (5.5 * c * ctx->p50 - 2.5 * ctx->p40) / 3.0;
            ctx->p11 = s;
            ctx->p21 = 3.0 * c * s;
            ctx->p31 = s * (7.5 * c2 - 1.5);
            ctx->p41 = c * s * (17.5 * c2 - 7.5);
            ctx->p51 = s * (39.375 * c4 - 26.25 * c2 + 1.875);
            ctx->p22 = 3.0 * s2;
            ctx->p32 = 15.0 * c * s2;
            ctx->p42 = s2 * (52.5 * c2 - 7.5);
            ctx->p52 = 3.0 * c * ctx->p42 - 2.0 * ctx->p32;
            ctx->p62 = 2.75 * c * ctx->p52 - 1.75 * ctx->p42;
            ctx->p33 = 15.0 * s * s2;

            // compute Legendre polynomials wrt magnetic pole (79N, 71W)
            double clmlmg = cos(lon - XLMG);
            double cmg  = s * CPMG * clmlmg + c * SPMG;
            double cmg2 = cmg * cmg;
            double cmg4 = cmg2 * cmg2;
            ctx->p10mg = cmg;
            ctx->p20mg = 1.5 * cmg2 - 0.5;
            ctx->p40mg = 4.375 * cmg4 - 3.75 * cmg2 + 0.375;

            ctx->clfl = scLon_c;
            ctx->slfl = scLon_s;

            // local time
            ctx->hl0 = hl;
            double scHlo_s=sin(ctx->hl0),scHlo_c=cos(ctx->hl0);
            ctx->ch  = scHlo_c;
            ctx->sh  = scHlo_s;
            ctx->c2h = ctx->ch * ctx->ch - ctx->sh * ctx->sh;
            ctx->s2h = 2.0 * ctx->ch * ctx->sh;
            ctx->c3h = ctx->c2h * ctx->ch - ctx->s2h * ctx->sh;
            ctx->s3h = ctx->s2h * ctx->ch + ctx->c2h * ctx->sh;

            double zlb = ZLB0; // + dzlb ??





















            //  compute function g(l) / tinf, t120, tp120
            int kleq = 1;
            double gdelt = g_function(ctx,dtm_tt, da, 1, kleq);
            da[1] = 1.0 + gdelt;
            double tinf   = dtm_tt[1] * da[1];

            kleq = 0; // equinox

            if (ctx->day < 59 || ctx->day > 284) {
                kleq = -1; // north winter
            }
            if (ctx->day > 99 && ctx->day < 244) {
                kleq = 1; // north summer
            }

            double gdelt0 =  g_function(ctx,dtm_t0, da, 0, kleq);
            da[1] = (dtm_t0[1] + gdelt0) / dtm_t0[1];
            double t120 = dtm_t0[1] + gdelt0;
            double gdeltp = g_function(ctx,dtm_tp, da, 0, kleq);
            da[1] = (dtm_tp[1] + gdeltp) / dtm_tp[1];
            double tp120 = dtm_tp[1] + gdeltp;

            // compute n(z) concentrations: H, He, O, N2, O2, N
            double sigma   = tp120 / (tinf - t120);
            double dzeta   = (RE + zlb) / (RE + altiKM);
            double zeta    = (altiKM - zlb) * dzeta;
            double sigzeta = sigma * zeta;
            double expsz   = exp(-sigzeta);
            double tz = tinf - (tinf - t120) * expsz;

            double dbase[7]={0};

            kleq = 1;

            double gdelh = g_function(ctx,dtm_h, da, 0, kleq);
            da[1] = exp(gdelh);
            dbase[1] = dtm_h[1] * da[1];

            double gdelhe = g_function(ctx,dtm_he, da, 0, kleq);
            da[1] = exp(gdelhe);
            dbase[2] = dtm_he[1] * da[1];

            double gdelo = g_function(ctx,dtm_o, da, 1, kleq);
            da[1] = exp(gdelo);
            dbase[3] = dtm_o[1] * da[1];

            double gdelaz2 = g_function(ctx,dtm_az2, da, 1, kleq);
            da[1] = exp(gdelaz2);
            dbase[4] = dtm_az2[1] * da[1];

            double gdelo2 = g_function(ctx,dtm_o2, da, 1, kleq);
            da[1] = exp(gdelo2);
            dbase[5] = dtm_o2[1] * da[1];

            double gdelaz = g_function(ctx,dtm_az, da, 1, kleq);
            da[1] = exp(gdelaz);
            dbase[6] = dtm_az[1] * da[1];

            double zlbre  = 1.0 + zlb / RE;
            double glb    = (GSURF / (zlbre * zlbre)) / (sigma * RGAS * tinf);
            double t120tz = t120 / tz;

            // Partial densities in (g/cm3).
            // d(1) = hydrogen
            // d(2) = helium
            // d(3) = atomic oxygen
            // d(4) = molecular nitrogen
            // d(5) = molecular oxygen
            // d(6) = atomic nitrogen
            double tmpro = 0;
            for (int i = 1; i <= 6; i++) {
                double gamma = MA[i] * glb;
                double upapg = 1.0 + ALEFA[i] + gamma;
                double fzI = exp(log(t120tz) * upapg - sigzeta * gamma);
                // concentrations of H, He, O, N2, O2, N (particles/cm³)
                double ccI = dbase[i] * fzI;
                // contribution of densities of H, He, O, N2, O2, N (g/cm³)
                tmpro += ccI * VMA[i];
                out->number_density_m3[i-1]=ccI*1e6;
                out->partial_mass_density_kg_m3[i-1]=ccI*VMA[i]*1000;
            }
            out->density_kg_m3=tmpro*1000;out->temperature_k=tz;out->exospheric_temperature_k=tinf;


}

tinydtm_status tinydtm_evaluate(const tinydtm_input *in,tinydtm_workspace *ctx,tinydtm_output *out) {
    tinydtm_output r;int i;double lon,hl;
    if(!in||!ctx||!out)return TINYDTM_INVALID;
    if(!isfinite(in->altitude_m)||!isfinite(in->longitude_rad)||!isfinite(in->latitude_rad)||!isfinite(in->local_solar_time_rad)||!isfinite(in->f107)||!isfinite(in->f107_mean)||!isfinite(in->kp_3h)||!isfinite(in->kp_24h))return TINYDTM_INVALID;
    if(in->day_of_year<1||in->day_of_year>366||in->altitude_m<120000||in->altitude_m>1500000||fabs(in->latitude_rad)>1.57079632679489661923||in->f107<40||in->f107>400||in->f107_mean<40||in->f107_mean>400||in->kp_3h<0||in->kp_3h>9||in->kp_24h<0||in->kp_24h>9)return TINYDTM_RANGE;
    ctx->day=in->day_of_year;ctx->f[0]=0;ctx->f[1]=in->f107;ctx->f[2]=0;ctx->fbar[0]=0;ctx->fbar[1]=in->f107_mean;ctx->fbar[2]=0;
    for(i=0;i<5;i++)ctx->akp[i]=0;
    ctx->akp[1]=in->kp_3h;ctx->akp[3]=in->kp_24h;
    for(i=0;i<97;i++)ctx->scratch[i]=0;
    lon=remainder(in->longitude_rad,6.28318530717958647693);hl=remainder(in->local_solar_time_rad,6.28318530717958647693);
    compute(ctx,in->altitude_m/1000,lon,in->latitude_rad,hl,&r);
    if(!isfinite(r.density_kg_m3)||r.density_kg_m3<0||!isfinite(r.temperature_k)||r.temperature_k<=0||!isfinite(r.exospheric_temperature_k)||r.exospheric_temperature_k<=0)return TINYDTM_ARITHMETIC;
    for(i=0;i<6;i++)if(!isfinite(r.number_density_m3[i])||r.number_density_m3[i]<0||!isfinite(r.partial_mass_density_kg_m3[i])||r.partial_mass_density_kg_m3[i]<0)return TINYDTM_ARITHMETIC;
    *out=r;return TINYDTM_OK;
}
