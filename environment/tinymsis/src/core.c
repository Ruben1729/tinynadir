/* TinyNadir changes: caller context, fixed spline workspace, immutable data,
   no diagnostic I/O, private names; original public-domain model below. */
/* -------------------------------------------------------------------- */
/* ---------  N R L M S I S E - 0 0    M O D E L    2 0 0 1  ---------- */
/* -------------------------------------------------------------------- */

/* This file is part of the NRLMSISE-00  C source code package - release
 * 20041227
 *
 * The NRLMSISE-00 model was developed by Mike Picone, Alan Hedin, and
 * Doug Drob. They also wrote a NRLMSISE-00 distribution package in
 * FORTRAN which is available at
 * http://uap-www.nrl.navy.mil/models_web/msis/msis_home.htm
 *
 * Dominik Brodowski implemented and maintains this C version. You can
 * reach him at mail@brodo.de. See the file "DOCUMENTATION" for details,
 * and check http://www.brodo.de/english/pub/nrlmsise/index.html for
 * updated releases of this package.
 */



/* ------------------------------------------------------------------- */
/* ------------------------------ INCLUDES --------------------------- */
/* ------------------------------------------------------------------- */

#include "core.h"
static void tinymsis_core_tselec(tinymsis_workspace *ctx, struct nrlmsise_flags *flags);
static void tinymsis_core_glatf(tinymsis_workspace *ctx, double lat, double *gv, double *reff);
static double tinymsis_core_ccor(tinymsis_workspace *ctx, double alt, double r, double h1, double zh);
static double tinymsis_core_ccor2(tinymsis_workspace *ctx, double alt, double r, double h1, double zh, double h2);
static double tinymsis_core_scalh(tinymsis_workspace *ctx, double alt, double xm, double temp);
static double tinymsis_core_dnet(tinymsis_workspace *ctx, double diffusive_density, double dm, double zhm, double xmm, double xm);
static void tinymsis_core_splini(tinymsis_workspace *ctx, double *xa, double *ya, double *y2a, int n, double x, double *y);
static void tinymsis_core_splint(tinymsis_workspace *ctx, double *xa, double *ya, double *y2a, int n, double x, double *y);
static void tinymsis_core_spline(tinymsis_workspace *ctx, double *x, double *y, int n, double yp1, double ypn, double *y2);
static double tinymsis_core_zeta(tinymsis_workspace *ctx, double zz, double zl);
static double tinymsis_core_densm(tinymsis_workspace *ctx, double alt, double d0, double xm, double *tz, int mn3, double *zn3, double *tn3, double *tgn3, int mn2, double *zn2, double *tn2, double *tgn2);
static double tinymsis_core_densu(tinymsis_workspace *ctx, double alt, double dlb, double tinf, double tlb, double xm, double alpha, double *tz, double zlb, double s2, int mn1, double *zn1, double *tn1, double *tgn1);
static double tinymsis_core_g0(tinymsis_workspace *ctx, double a, const double *p);
static double tinymsis_core_sumex(tinymsis_workspace *ctx, double ex);
static double tinymsis_core_sg0(tinymsis_workspace *ctx, double ex, const double *p, double *ap);
static double tinymsis_core_globe7(tinymsis_workspace *ctx, const double *p, struct nrlmsise_input *input, struct nrlmsise_flags *flags);
static double tinymsis_core_glob7s(tinymsis_workspace *ctx, const double *p, struct nrlmsise_input *input, struct nrlmsise_flags *flags);
static void tinymsis_core_gts7(tinymsis_workspace *ctx, struct nrlmsise_input *input, struct nrlmsise_flags *flags, struct nrlmsise_output *output);   /* header for nrlmsise-00.h */
#include <math.h>          /* maths functions */





/* ------------------------------------------------------------------- */
/* ------------------------- SHARED VARIABLES ------------------------ */
/* ------------------------------------------------------------------- */

/* PARMB */

/* GTS3C */

/* DMIX */

/* MESO7 */

/* POWER7 */
extern const double tinymsis_coeff_pt[150];
extern const double tinymsis_coeff_pd[9][150];
extern const double tinymsis_coeff_ps[150];
extern const double tinymsis_coeff_pdl[2][25];
extern const double tinymsis_coeff_ptl[4][100];
extern const double tinymsis_coeff_pma[10][100];
extern const double tinymsis_coeff_sam[100];

/* LOWER7 */
extern const double tinymsis_coeff_ptm[10];
extern const double tinymsis_coeff_pdm[8][10];
extern const double tinymsis_coeff_pavgm[10];

/* LPOLY */



/* ------------------------------------------------------------------- */
/* ------------------------------ TSELEC ----------------------------- */
/* ------------------------------------------------------------------- */

static void tinymsis_core_tselec(tinymsis_workspace *ctx, struct nrlmsise_flags *flags) {
    (void)ctx;
    int i;
    for (i=0;i<24;i++) {
        if (i!=9) {
            if (flags->switches[i]==1)
                flags->sw[i]=1;
            else
                flags->sw[i]=0;
            if (flags->switches[i]>0)
                flags->swc[i]=1;
            else
                flags->swc[i]=0;
        } else {
            flags->sw[i]=flags->switches[i];
            flags->swc[i]=flags->switches[i];
        }
    }
}



/* ------------------------------------------------------------------- */
/* ------------------------------ GLATF ------------------------------ */
/* ------------------------------------------------------------------- */

static void tinymsis_core_glatf(tinymsis_workspace *ctx, double lat, double *gv, double *reff) {
    (void)ctx;
    double dgtr = 1.74533E-2;
    double c2;
    c2 = cos(2.0*dgtr*lat);
    *gv = 980.616 * (1.0 - 0.0026373 * c2);
    *reff = 2.0 * (*gv) / (3.085462E-6 + 2.27E-9 * c2) * 1.0E-5;
}



/* ------------------------------------------------------------------- */
/* ------------------------------ CCOR ------------------------------- */
/* ------------------------------------------------------------------- */

static double tinymsis_core_ccor(tinymsis_workspace *ctx, double alt, double r, double h1, double zh) {
    (void)ctx;
/*        CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
 *         ALT - altitude
 *         R - target ratio
 *         H1 - transition scale length
 *         ZH - altitude of 1/2 R
 */
    double e;
    double ex;
    e = (alt - zh) / h1;
    if (e>70)
        return exp(0);
    if (e<-70)
        return exp(r);
    ex = exp(e);
    e = r / (1.0 + ex);
    return exp(e);
}



/* ------------------------------------------------------------------- */
/* ------------------------------ CCOR ------------------------------- */
/* ------------------------------------------------------------------- */

static double tinymsis_core_ccor2(tinymsis_workspace *ctx, double alt, double r, double h1, double zh, double h2) {
    (void)ctx;
/*        CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
 *         ALT - altitude
 *         R - target ratio
 *         H1 - transition scale length
 *         ZH - altitude of 1/2 R
 *         H2 - transition scale length #2 ?
 */
    double e1, e2;
    double ex1, ex2;
    double ccor2v;
    e1 = (alt - zh) / h1;
    e2 = (alt - zh) / h2;
    if ((e1 > 70) || (e2 > 70))
        return exp(0);
    if ((e1 < -70) && (e2 < -70))
        return exp(r);
    ex1 = exp(e1);
    ex2 = exp(e2);
    ccor2v = r / (1.0 + 0.5 * (ex1 + ex2));
    return exp(ccor2v);
}



/* ------------------------------------------------------------------- */
/* ------------------------------- SCALH ----------------------------- */
/* ------------------------------------------------------------------- */

static double tinymsis_core_scalh(tinymsis_workspace *ctx, double alt, double xm, double temp) {
    (void)ctx;
    double g;
    double rgas=831.4;
    g = ctx->gsurf / (pow((1.0 + alt/ctx->re),2.0));
    g = rgas * temp / (g * xm);
    return g;
}



/* ------------------------------------------------------------------- */
/* -------------------------------- DNET ----------------------------- */
/* ------------------------------------------------------------------- */

static double tinymsis_core_dnet(tinymsis_workspace *ctx, double diffusive_density, double dm, double zhm, double xmm, double xm) {
    (void)ctx;
/*       TURBOPAUSE CORRECTION FOR MSIS MODELS
 *        Root mean density
 *         DD - diffusive density
 *         DM - full mixed density
 *         ZHM - transition scale length
 *         XMM - full mixed molecular weight
 *         XM  - species molecular weight
 *         DNET - combined density
 */
    double a;
    double ylog;
    a  = zhm / (xmm-xm);
    if (!((dm>0) && (diffusive_density>0))) {

        if ((diffusive_density==0) && (dm==0))
            diffusive_density=1;
        if (dm==0)
            return diffusive_density;
        if (diffusive_density==0)
            return dm;
    }
    ylog = a * log(dm/diffusive_density);
    if (ylog<-10)
        return diffusive_density;
    if (ylog>10)
        return dm;
    a = diffusive_density*pow((1.0 + exp(ylog)),(1.0/a));
    return a;
}



/* ------------------------------------------------------------------- */
/* ------------------------------- SPLINI ---------------------------- */
/* ------------------------------------------------------------------- */

static void tinymsis_core_splini(tinymsis_workspace *ctx, double *xa, double *ya, double *y2a, int n, double x, double *y) {
    (void)ctx;
/*      INTEGRATE CUBIC SPLINE FUNCTION FROM XA(1) TO X
 *       XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
 *       Y2A: ARRAY OF SECOND DERIVATIVES
 *       N: SIZE OF ARRAYS XA,YA,Y2A
 *       X: ABSCISSA ENDPOINT FOR INTEGRATION
 *       Y: OUTPUT VALUE
 */
    double yi=0;
    int klo=0;
    int khi=1;
    double xx, h, a, b, a2, b2;
    while ((x>xa[klo]) && (khi<n)) {
        xx=x;
        if (khi<(n-1)) {
            if (x<xa[khi])
                xx=x;
            else
                xx=xa[khi];
        }
        h = xa[khi] - xa[klo];
        a = (xa[khi] - xx)/h;
        b = (xx - xa[klo])/h;
        a2 = a*a;
        b2 = b*b;
        yi += ((1.0 - a2) * ya[klo] / 2.0 + b2 * ya[khi] / 2.0 + ((-(1.0+a2*a2)/4.0 + a2/2.0) * y2a[klo] + (b2*b2/4.0 - b2/2.0) * y2a[khi]) * h * h / 6.0) * h;
        klo++;
        khi++;
    }
    *y = yi;
}



/* ------------------------------------------------------------------- */
/* ------------------------------- SPLINT ---------------------------- */
/* ------------------------------------------------------------------- */

static void tinymsis_core_splint(tinymsis_workspace *ctx, double *xa, double *ya, double *y2a, int n, double x, double *y) {
    (void)ctx;
/*      CALCULATE CUBIC SPLINE INTERP VALUE
 *       ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL.
 *       XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
 *       Y2A: ARRAY OF SECOND DERIVATIVES
 *       N: SIZE OF ARRAYS XA,YA,Y2A
 *       X: ABSCISSA FOR INTERPOLATION
 *       Y: OUTPUT VALUE
 */
    int klo=0;
    int khi=n-1;
    int k;
    double h;
    double a, b, yi;
    while ((khi-klo)>1) {
        k=(khi+klo)/2;
        if (xa[k]>x)
            khi=k;
        else
            klo=k;
    }
    h = xa[khi] - xa[klo];
    if (h==0.0)
        ctx->numeric_error=1;
    a = (xa[khi] - x)/h;
    b = (x - xa[klo])/h;
    yi = a * ya[klo] + b * ya[khi] + ((a*a*a - a) * y2a[klo] + (b*b*b - b) * y2a[khi]) * h * h/6.0;
    *y = yi;
}



/* ------------------------------------------------------------------- */
/* ------------------------------- SPLINE ---------------------------- */
/* ------------------------------------------------------------------- */

static void tinymsis_core_spline(tinymsis_workspace *ctx, double *x, double *y, int n, double yp1, double ypn, double *y2) {
    (void)ctx;
/*       CALCULATE 2ND DERIVATIVES OF CUBIC SPLINE INTERP FUNCTION
 *       ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL
 *       X,Y: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
 *       N: SIZE OF ARRAYS X,Y
 *       YP1,YPN: SPECIFIED DERIVATIVES AT X[0] AND X[N-1]; VALUES
 *                >= 1E30 SIGNAL SIGNAL SECOND DERIVATIVE ZERO
 *       Y2: OUTPUT ARRAY OF SECOND DERIVATIVES
 */
    double u[10];
    double sig, p, qn, un;
    int i, k;
    if(n<2||n>10){ctx->numeric_error=1;return;}
    if (yp1>0.99E30) {
        y2[0]=0;
        u[0]=0;
    } else {
        y2[0]=-0.5;
        u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
    }
    for (i=1;i<(n-1);i++) {
        sig = (x[i]-x[i-1])/(x[i+1] - x[i-1]);
        p = sig * y2[i-1] + 2.0;
        y2[i] = (sig - 1.0) / p;
        u[i] = (6.0 * ((y[i+1] - y[i])/(x[i+1] - x[i]) -(y[i] - y[i-1]) / (x[i] - x[i-1]))/(x[i+1] - x[i-1]) - sig * u[i-1])/p;
    }
    if (ypn>0.99E30) {
        qn = 0;
        un = 0;
    } else {
        qn = 0.5;
        un = (3.0 / (x[n-1] - x[n-2])) * (ypn - (y[n-1] - y[n-2])/(x[n-1] - x[n-2]));
    }
    y2[n-1] = (un - qn * u[n-2]) / (qn * y2[n-2] + 1.0);
    for (k=n-2;k>=0;k--)
        y2[k] = y2[k] * y2[k+1] + u[k];


}



/* ------------------------------------------------------------------- */
/* ------------------------------- DENSM ----------------------------- */
/* ------------------------------------------------------------------- */

static double tinymsis_core_zeta(tinymsis_workspace *ctx, double zz, double zl) {
    (void)ctx;
    return ((zz-zl)*(ctx->re+zl)/(ctx->re+zz));
}

static double tinymsis_core_densm(tinymsis_workspace *ctx, double alt, double d0, double xm, double *tz, int mn3, double *zn3, double *tn3, double *tgn3, int mn2, double *zn2, double *tn2, double *tgn2) {
    (void)ctx;
/*      Calculate Temperature and Density Profiles for lower atmos.  */
    double xs[10], ys[10], y2out[10];
    double rgas = 831.4;
    double z, z1, z2, t1, t2, zg, zgdif;
    double yd1, yd2;
    double x, y, yi;
    double expl, gamm, glb;
    double densm_tmp;
    int mn;
    int k;
    densm_tmp=d0;
    if (alt>zn2[0]) {
        if (xm==0.0)
            return *tz;
        else
            return d0;
    }

    /* STRATOSPHERE/MESOSPHERE TEMPERATURE */
    if (alt>zn2[mn2-1])
        z=alt;
    else
        z=zn2[mn2-1];
    mn=mn2;
    z1=zn2[0];
    z2=zn2[mn-1];
    t1=tn2[0];
    t2=tn2[mn-1];
    zg = tinymsis_core_zeta(ctx, z, z1);
    zgdif = tinymsis_core_zeta(ctx, z2, z1);

    /* set up tinymsis_core_spline nodes */
    for (k=0;k<mn;k++) {
        xs[k]=tinymsis_core_zeta(ctx, zn2[k],z1)/zgdif;
        ys[k]=1.0 / tn2[k];
    }
    yd1=-tgn2[0] / (t1*t1) * zgdif;
    yd2=-tgn2[1] / (t2*t2) * zgdif * (pow(((ctx->re+z2)/(ctx->re+z1)),2.0));

    /* calculate tinymsis_core_spline coefficients */
    tinymsis_core_spline(ctx, xs, ys, mn, yd1, yd2, y2out);
    x = zg/zgdif;
    tinymsis_core_splint(ctx, xs, ys, y2out, mn, x, &y);

    /* temperature at altitude */
    *tz = 1.0 / y;
    if (xm!=0.0) {
        /* calaculate stratosphere / mesospehere density */
        glb = ctx->gsurf / (pow((1.0 + z1/ctx->re),2.0));
        gamm = xm * glb * zgdif / rgas;

        /* Integrate temperature profile */
        tinymsis_core_splini(ctx, xs, ys, y2out, mn, x, &yi);
        expl=gamm*yi;
        if (expl>50.0)
            expl=50.0;

        /* Density at altitude */
        densm_tmp = densm_tmp * (t1 / *tz) * exp(-expl);
    }

    if (alt>zn3[0]) {
        if (xm==0.0)
            return *tz;
        else
            return densm_tmp;
    }

    /* troposhere / stratosphere temperature */
    z = alt;
    mn = mn3;
    z1=zn3[0];
    z2=zn3[mn-1];
    t1=tn3[0];
    t2=tn3[mn-1];
    zg=tinymsis_core_zeta(ctx, z,z1);
    zgdif=tinymsis_core_zeta(ctx, z2,z1);

    /* set up tinymsis_core_spline nodes */
    for (k=0;k<mn;k++) {
        xs[k] = tinymsis_core_zeta(ctx, zn3[k],z1) / zgdif;
        ys[k] = 1.0 / tn3[k];
    }
    yd1=-tgn3[0] / (t1*t1) * zgdif;
    yd2=-tgn3[1] / (t2*t2) * zgdif * (pow(((ctx->re+z2)/(ctx->re+z1)),2.0));

    /* calculate tinymsis_core_spline coefficients */
    tinymsis_core_spline(ctx, xs, ys, mn, yd1, yd2, y2out);
    x = zg/zgdif;
    tinymsis_core_splint(ctx, xs, ys, y2out, mn, x, &y);

    /* temperature at altitude */
    *tz = 1.0 / y;
    if (xm!=0.0) {
        /* calaculate tropospheric / stratosphere density */
        glb = ctx->gsurf / (pow((1.0 + z1/ctx->re),2.0));
        gamm = xm * glb * zgdif / rgas;

        /* Integrate temperature profile */
        tinymsis_core_splini(ctx, xs, ys, y2out, mn, x, &yi);
        expl=gamm*yi;
        if (expl>50.0)
            expl=50.0;

        /* Density at altitude */
        densm_tmp = densm_tmp * (t1 / *tz) * exp(-expl);
    }
    if (xm==0.0)
        return *tz;
    else
        return densm_tmp;
}



/* ------------------------------------------------------------------- */
/* ------------------------------- DENSU ----------------------------- */
/* ------------------------------------------------------------------- */

static double tinymsis_core_densu(tinymsis_workspace *ctx, double alt, double dlb, double tinf, double tlb, double xm, double alpha, double *tz, double zlb, double s2, int mn1, double *zn1, double *tn1, double *tgn1) {
    (void)ctx;
/*      Calculate Temperature and Density Profiles for MSIS models
 *      New lower thermo polynomial
 */
    double yd2, yd1, x=0, y;
    double rgas=831.4;
    double densu_temp=1.0;
    double za, z, zg2, tt, ta;
    double dta, z1=0, z2, t1=0, t2, zg, zgdif=0;
    int mn=0;
    int k;
    double glb;
    double expl;
    double yi;
    double densa;
    double gamma, gamm;
    double xs[5], ys[5], y2out[5];
    /* joining altitudes of Bates and tinymsis_core_spline */
    za=zn1[0];
    if (alt>za)
        z=alt;
    else
        z=za;

    /* geopotential altitude difference from ZLB */
    zg2 = tinymsis_core_zeta(ctx, z, zlb);

    /* Bates temperature */
    tt = tinf - (tinf - tlb) * exp(-s2*zg2);
    ta = tt;
    *tz = tt;
    densu_temp = *tz;

    if (alt<za) {
        /* calculate temperature below ZA
         * temperature gradient at ZA from Bates profile */
        dta = (tinf - ta) * s2 * pow(((ctx->re+zlb)/(ctx->re+za)),2.0);
        tgn1[0]=dta;
        tn1[0]=ta;
        if (alt>zn1[mn1-1])
            z=alt;
        else
            z=zn1[mn1-1];
        mn=mn1;
        z1=zn1[0];
        z2=zn1[mn-1];
        t1=tn1[0];
        t2=tn1[mn-1];
        /* geopotental difference from z1 */
        zg = tinymsis_core_zeta(ctx, z, z1);
        zgdif = tinymsis_core_zeta(ctx, z2, z1);
        /* set up tinymsis_core_spline nodes */
        for (k=0;k<mn;k++) {
            xs[k] = tinymsis_core_zeta(ctx, zn1[k], z1) / zgdif;
            ys[k] = 1.0 / tn1[k];
        }
        /* end node derivatives */
        yd1 = -tgn1[0] / (t1*t1) * zgdif;
        yd2 = -tgn1[1] / (t2*t2) * zgdif * pow(((ctx->re+z2)/(ctx->re+z1)),2.0);
        /* calculate tinymsis_core_spline coefficients */
        tinymsis_core_spline(ctx, xs, ys, mn, yd1, yd2, y2out);
        x = zg / zgdif;
        tinymsis_core_splint(ctx, xs, ys, y2out, mn, x, &y);
        /* temperature at altitude */
        *tz = 1.0 / y;
        densu_temp = *tz;
    }
    if (xm==0)
        return densu_temp;

    /* calculate density above za */
    glb = ctx->gsurf / pow((1.0 + zlb/ctx->re),2.0);
    gamma = xm * glb / (s2 * rgas * tinf);
    expl = exp(-s2 * gamma * zg2);
    if (expl>50.0)
        expl=50.0;
    if (tt<=0)
        expl=50.0;

    /* density at altitude */
    densa = dlb * pow((tlb/tt),((1.0+alpha+gamma))) * expl;
    densu_temp=densa;
    if (alt>=za)
        return densu_temp;

    /* calculate density below za */
    glb = ctx->gsurf / pow((1.0 + z1/ctx->re),2.0);
    gamm = xm * glb * zgdif / rgas;

    /* integrate tinymsis_core_spline temperatures */
    tinymsis_core_splini(ctx, xs, ys, y2out, mn, x, &yi);
    expl = gamm * yi;
    if (expl>50.0)
        expl=50.0;
    if (*tz<=0)
        expl=50.0;

    /* density at altitude */
    densu_temp = densu_temp * pow ((t1 / *tz),(1.0 + alpha)) * exp(-expl);
    return densu_temp;
}



/* ------------------------------------------------------------------- */
/* ------------------------------- GLOBE7 ---------------------------- */
/* ------------------------------------------------------------------- */

/*    3hr Magnetic activity functions */
/*    Eq. A24d */
static double tinymsis_core_g0(tinymsis_workspace *ctx, double a, const double *p) {
    (void)ctx;
    double local_p24=fmax(p[24],1.0E-4);
    return (a - 4.0 + (p[25] - 1.0) * (a - 4.0 + (exp(-sqrt(local_p24*local_p24) * (a - 4.0)) - 1.0) / sqrt(local_p24*local_p24)));
}

/*    Eq. A24c */
static double tinymsis_core_sumex(tinymsis_workspace *ctx, double ex) {
    (void)ctx;
    return (1.0 + (1.0 - pow(ex,19.0)) / (1.0 - ex) * pow(ex,0.5));
}

/*    Eq. A24a */
static double tinymsis_core_sg0(tinymsis_workspace *ctx, double ex, const double *p, double *ap) {
    (void)ctx;
    return (tinymsis_core_g0(ctx, ap[1],p) + (tinymsis_core_g0(ctx, ap[2],p)*ex + tinymsis_core_g0(ctx, ap[3],p)*ex*ex + \
                tinymsis_core_g0(ctx, ap[4],p)*pow(ex,3.0)	+ (tinymsis_core_g0(ctx, ap[5],p)*pow(ex,4.0) + \
                tinymsis_core_g0(ctx, ap[6],p)*pow(ex,12.0))*(1.0-pow(ex,8.0))/(1.0-ex)))/tinymsis_core_sumex(ctx, ex);
}

static double tinymsis_core_globe7(tinymsis_workspace *ctx, const double *p, struct nrlmsise_input *input, struct nrlmsise_flags *flags) {
    (void)ctx;
/*       CALCULATE G(L) FUNCTION
 *       Upper Thermosphere Parameters */
    double t[15];
    int i,j;
    double apd;
    double tloc;
    double c, s, c2, c4, s2;
    double sr = 7.2722E-5;
    double dgtr = 1.74533E-2;
    double dr = 1.72142E-2;
    double hr = 0.2618;
    double cd32, cd18, cd14, cd39;
    double df;
    double f1, f2;
    double tinf;
    struct ap_array *ap;

    tloc=input->lst;
    for (j=0;j<14;j++)
        t[j]=0;

    /* calculate legendre polynomials */
    c = sin(input->g_lat * dgtr);
    s = cos(input->g_lat * dgtr);
    c2 = c*c;
    c4 = c2*c2;
    s2 = s*s;

    ctx->plg[0][1] = c;
    ctx->plg[0][2] = 0.5*(3.0*c2 -1.0);
    ctx->plg[0][3] = 0.5*(5.0*c*c2-3.0*c);
    ctx->plg[0][4] = (35.0*c4 - 30.0*c2 + 3.0)/8.0;
    ctx->plg[0][5] = (63.0*c2*c2*c - 70.0*c2*c + 15.0*c)/8.0;
    ctx->plg[0][6] = (11.0*c*ctx->plg[0][5] - 5.0*ctx->plg[0][4])/6.0;
/*      ctx->plg[0][7] = (13.0*c*ctx->plg[0][6] - 6.0*ctx->plg[0][5])/7.0; */
    ctx->plg[1][1] = s;
    ctx->plg[1][2] = 3.0*c*s;
    ctx->plg[1][3] = 1.5*(5.0*c2-1.0)*s;
    ctx->plg[1][4] = 2.5*(7.0*c2*c-3.0*c)*s;
    ctx->plg[1][5] = 1.875*(21.0*c4 - 14.0*c2 +1.0)*s;
    ctx->plg[1][6] = (11.0*c*ctx->plg[1][5]-6.0*ctx->plg[1][4])/5.0;
/*      ctx->plg[1][7] = (13.0*c*ctx->plg[1][6]-7.0*ctx->plg[1][5])/6.0; */
/*      ctx->plg[1][8] = (15.0*c*ctx->plg[1][7]-8.0*ctx->plg[1][6])/7.0; */
    ctx->plg[2][2] = 3.0*s2;
    ctx->plg[2][3] = 15.0*s2*c;
    ctx->plg[2][4] = 7.5*(7.0*c2 -1.0)*s2;
    ctx->plg[2][5] = 3.0*c*ctx->plg[2][4]-2.0*ctx->plg[2][3];
    ctx->plg[2][6] =(11.0*c*ctx->plg[2][5]-7.0*ctx->plg[2][4])/4.0;
    ctx->plg[2][7] =(13.0*c*ctx->plg[2][6]-8.0*ctx->plg[2][5])/5.0;
    ctx->plg[3][3] = 15.0*s2*s;
    ctx->plg[3][4] = 105.0*s2*s*c;
    ctx->plg[3][5] =(9.0*c*ctx->plg[3][4]-7.*ctx->plg[3][3])/2.0;
    ctx->plg[3][6] =(11.0*c*ctx->plg[3][5]-8.*ctx->plg[3][4])/3.0;

    if (!(((flags->sw[7]==0)&&(flags->sw[8]==0))&&(flags->sw[14]==0))) {
        ctx->stloc = sin(hr*tloc);
        ctx->ctloc = cos(hr*tloc);
        ctx->s2tloc = sin(2.0*hr*tloc);
        ctx->c2tloc = cos(2.0*hr*tloc);
        ctx->s3tloc = sin(3.0*hr*tloc);
        ctx->c3tloc = cos(3.0*hr*tloc);
    }

    cd32 = cos(dr*(input->doy-p[31]));
    cd18 = cos(2.0*dr*(input->doy-p[17]));
    cd14 = cos(dr*(input->doy-p[13]));
    cd39 = cos(2.0*dr*(input->doy-p[38]));

    /* F10.7 EFFECT */
    df = input->f107 - input->f107A;
    ctx->dfa = input->f107A - 150.0;
    t[0] =  p[19]*df*(1.0+p[59]*ctx->dfa) + p[20]*df*df + p[21]*ctx->dfa + p[29]*pow(ctx->dfa,2.0);
    f1 = 1.0 + (p[47]*ctx->dfa +p[19]*df+p[20]*df*df)*flags->swc[1];
    f2 = 1.0 + (p[49]*ctx->dfa+p[19]*df+p[20]*df*df)*flags->swc[1];

    /*  TIME INDEPENDENT */
    t[1] = (p[1]*ctx->plg[0][2]+ p[2]*ctx->plg[0][4]+p[22]*ctx->plg[0][6]) + \
          (p[14]*ctx->plg[0][2])*ctx->dfa*flags->swc[1] +p[26]*ctx->plg[0][1];

    /*  SYMMETRICAL ANNUAL */
    t[2] = p[18]*cd32;

    /*  SYMMETRICAL SEMIANNUAL */
    t[3] = (p[15]+p[16]*ctx->plg[0][2])*cd18;

    /*  ASYMMETRICAL ANNUAL */
    t[4] =  f1*(p[9]*ctx->plg[0][1]+p[10]*ctx->plg[0][3])*cd14;

    /*  ASYMMETRICAL SEMIANNUAL */
    t[5] =    p[37]*ctx->plg[0][1]*cd39;

        /* DIURNAL */
    if (flags->sw[7]) {
        double t71, t72;
        t71 = (p[11]*ctx->plg[1][2])*cd14*flags->swc[5];
        t72 = (p[12]*ctx->plg[1][2])*cd14*flags->swc[5];
        t[6] = f2*((p[3]*ctx->plg[1][1] + p[4]*ctx->plg[1][3] + p[27]*ctx->plg[1][5] + t71) * \
               ctx->ctloc + (p[6]*ctx->plg[1][1] + p[7]*ctx->plg[1][3] + p[28]*ctx->plg[1][5] \
                    + t72)*ctx->stloc);
}

    /* SEMIDIURNAL */
    if (flags->sw[8]) {
        double t81, t82;
        t81 = (p[23]*ctx->plg[2][3]+p[35]*ctx->plg[2][5])*cd14*flags->swc[5];
        t82 = (p[33]*ctx->plg[2][3]+p[36]*ctx->plg[2][5])*cd14*flags->swc[5];
        t[7] = f2*((p[5]*ctx->plg[2][2]+ p[41]*ctx->plg[2][4] + t81)*ctx->c2tloc +(p[8]*ctx->plg[2][2] + p[42]*ctx->plg[2][4] + t82)*ctx->s2tloc);
    }

    /* TERDIURNAL */
    if (flags->sw[14]) {
        t[13] = f2 * ((p[39]*ctx->plg[3][3]+(p[93]*ctx->plg[3][4]+p[46]*ctx->plg[3][6])*cd14*flags->swc[5])* ctx->s3tloc +(p[40]*ctx->plg[3][3]+(p[94]*ctx->plg[3][4]+p[48]*ctx->plg[3][6])*cd14*flags->swc[5])* ctx->c3tloc);
}

    /* magnetic activity based on daily ap */
    if (flags->sw[9]==-1) {
        ap = input->ap_a;
        if (p[51]!=0) {
            double exp1;
            exp1 = exp(-10800.0*sqrt(p[51]*p[51])/(1.0+p[138]*(45.0-sqrt(input->g_lat*input->g_lat))));
            if (exp1>0.99999)
                exp1=0.99999;
            ctx->apt[0]=tinymsis_core_sg0(ctx, exp1,p,ap->a);
            /* ctx->apt[1]=sg2(exp1,p,ap->a);
               ctx->apt[2]=tinymsis_core_sg0(ctx, exp2,p,ap->a);
               ctx->apt[3]=sg2(exp2,p,ap->a);
            */
            if (flags->sw[9]) {
                t[8] = ctx->apt[0]*(p[50]+p[96]*ctx->plg[0][2]+p[54]*ctx->plg[0][4]+ \
     (p[125]*ctx->plg[0][1]+p[126]*ctx->plg[0][3]+p[127]*ctx->plg[0][5])*cd14*flags->swc[5]+ \
     (p[128]*ctx->plg[1][1]+p[129]*ctx->plg[1][3]+p[130]*ctx->plg[1][5])*flags->swc[7]* \
                           cos(hr*(tloc-p[131])));
            }
        }
    } else {
        double p44, p45;
        apd=input->ap-4.0;
        p44=p[43];
        p45=p[44];
        if (p44<0)
            p44 = 1.0E-5;
        ctx->apdf = apd + (p45-1.0)*(apd + (exp(-p44 * apd) - 1.0)/p44);
        if (flags->sw[9]) {
            t[8]=ctx->apdf*(p[32]+p[45]*ctx->plg[0][2]+p[34]*ctx->plg[0][4]+ \
     (p[100]*ctx->plg[0][1]+p[101]*ctx->plg[0][3]+p[102]*ctx->plg[0][5])*cd14*flags->swc[5]+
     (p[121]*ctx->plg[1][1]+p[122]*ctx->plg[1][3]+p[123]*ctx->plg[1][5])*flags->swc[7]*
                    cos(hr*(tloc-p[124])));
        }
    }

    if ((flags->sw[10])&&(input->g_long>-1000.0)) {

        /* longitudinal */
        if (flags->sw[11]) {
            t[10] = (1.0 + p[80]*ctx->dfa*flags->swc[1])* \
     ((p[64]*ctx->plg[1][2]+p[65]*ctx->plg[1][4]+p[66]*ctx->plg[1][6]\
      +p[103]*ctx->plg[1][1]+p[104]*ctx->plg[1][3]+p[105]*ctx->plg[1][5]\
      +flags->swc[5]*(p[109]*ctx->plg[1][1]+p[110]*ctx->plg[1][3]+p[111]*ctx->plg[1][5])*cd14)* \
          cos(dgtr*input->g_long) \
      +(p[90]*ctx->plg[1][2]+p[91]*ctx->plg[1][4]+p[92]*ctx->plg[1][6]\
      +p[106]*ctx->plg[1][1]+p[107]*ctx->plg[1][3]+p[108]*ctx->plg[1][5]\
      +flags->swc[5]*(p[112]*ctx->plg[1][1]+p[113]*ctx->plg[1][3]+p[114]*ctx->plg[1][5])*cd14)* \
      sin(dgtr*input->g_long));
        }

        /* ut and mixed ut, longitude */
        if (flags->sw[12]){
            t[11]=(1.0+p[95]*ctx->plg[0][1])*(1.0+p[81]*ctx->dfa*flags->swc[1])*\
                (1.0+p[119]*ctx->plg[0][1]*flags->swc[5]*cd14)*\
                ((p[68]*ctx->plg[0][1]+p[69]*ctx->plg[0][3]+p[70]*ctx->plg[0][5])*\
                cos(sr*(input->sec-p[71])));
            t[11]+=flags->swc[11]*\
                (p[76]*ctx->plg[2][3]+p[77]*ctx->plg[2][5]+p[78]*ctx->plg[2][7])*\
                cos(sr*(input->sec-p[79])+2.0*dgtr*input->g_long)*(1.0+p[137]*ctx->dfa*flags->swc[1]);
        }

        /* ut, longitude magnetic activity */
        if (flags->sw[13]) {
            if (flags->sw[9]==-1) {
                if (p[51]) {
                    t[12]=ctx->apt[0]*flags->swc[11]*(1.+p[132]*ctx->plg[0][1])*\
                        ((p[52]*ctx->plg[1][2]+p[98]*ctx->plg[1][4]+p[67]*ctx->plg[1][6])*\
                         cos(dgtr*(input->g_long-p[97])))\
                        +ctx->apt[0]*flags->swc[11]*flags->swc[5]*\
                        (p[133]*ctx->plg[1][1]+p[134]*ctx->plg[1][3]+p[135]*ctx->plg[1][5])*\
                        cd14*cos(dgtr*(input->g_long-p[136])) \
                        +ctx->apt[0]*flags->swc[12]* \
                        (p[55]*ctx->plg[0][1]+p[56]*ctx->plg[0][3]+p[57]*ctx->plg[0][5])*\
                        cos(sr*(input->sec-p[58]));
                }
            } else {
                t[12] = ctx->apdf*flags->swc[11]*(1.0+p[120]*ctx->plg[0][1])*\
                    ((p[60]*ctx->plg[1][2]+p[61]*ctx->plg[1][4]+p[62]*ctx->plg[1][6])*\
                    cos(dgtr*(input->g_long-p[63])))\
                    +ctx->apdf*flags->swc[11]*flags->swc[5]* \
                    (p[115]*ctx->plg[1][1]+p[116]*ctx->plg[1][3]+p[117]*ctx->plg[1][5])* \
                    cd14*cos(dgtr*(input->g_long-p[118])) \
                    + ctx->apdf*flags->swc[12]* \
                    (p[83]*ctx->plg[0][1]+p[84]*ctx->plg[0][3]+p[85]*ctx->plg[0][5])* \
                    cos(sr*(input->sec-p[75]));
            }
        }
    }

    /* parms not used: 82, 89, 99, 139-149 */
    tinf = p[30];
    for (i=0;i<14;i++)
        tinf = tinf + fabs(flags->sw[i+1])*t[i];
    return tinf;
}



/* ------------------------------------------------------------------- */
/* ------------------------------- GLOB7S ---------------------------- */
/* ------------------------------------------------------------------- */

static double tinymsis_core_glob7s(tinymsis_workspace *ctx, const double *p, struct nrlmsise_input *input, struct nrlmsise_flags *flags) {
    (void)ctx;
/*    VERSION OF GLOBE FOR LOWER ATMOSPHERE 10/26/99
 */
    double pset=2.0;
    double t[14];
    double tt;
    double cd32, cd18, cd14, cd39;
    int i,j;
    double dr=1.72142E-2;
    double dgtr=1.74533E-2;
    /* confirm parameter set */
    if (p[99]!=pset && p[99]!=0) {
        ctx->numeric_error=1;
        return -1;
    }
    for (j=0;j<14;j++)
        t[j]=0.0;
    cd32 = cos(dr*(input->doy-p[31]));
    cd18 = cos(2.0*dr*(input->doy-p[17]));
    cd14 = cos(dr*(input->doy-p[13]));
    cd39 = cos(2.0*dr*(input->doy-p[38]));

    /* F10.7 */
    t[0] = p[21]*ctx->dfa;

    /* time independent */
    t[1]=p[1]*ctx->plg[0][2] + p[2]*ctx->plg[0][4] + p[22]*ctx->plg[0][6] + p[26]*ctx->plg[0][1] + p[14]*ctx->plg[0][3] + p[59]*ctx->plg[0][5];

        /* SYMMETRICAL ANNUAL */
    t[2]=(p[18]+p[47]*ctx->plg[0][2]+p[29]*ctx->plg[0][4])*cd32;

        /* SYMMETRICAL SEMIANNUAL */
    t[3]=(p[15]+p[16]*ctx->plg[0][2]+p[30]*ctx->plg[0][4])*cd18;

        /* ASYMMETRICAL ANNUAL */
    t[4]=(p[9]*ctx->plg[0][1]+p[10]*ctx->plg[0][3]+p[20]*ctx->plg[0][5])*cd14;

    /* ASYMMETRICAL SEMIANNUAL */
    t[5]=(p[37]*ctx->plg[0][1])*cd39;

        /* DIURNAL */
    if (flags->sw[7]) {
        double t71, t72;
        t71 = p[11]*ctx->plg[1][2]*cd14*flags->swc[5];
        t72 = p[12]*ctx->plg[1][2]*cd14*flags->swc[5];
        t[6] = ((p[3]*ctx->plg[1][1] + p[4]*ctx->plg[1][3] + t71) * ctx->ctloc + (p[6]*ctx->plg[1][1] + p[7]*ctx->plg[1][3] + t72) * ctx->stloc) ;
    }

    /* SEMIDIURNAL */
    if (flags->sw[8]) {
        double t81, t82;
        t81 = (p[23]*ctx->plg[2][3]+p[35]*ctx->plg[2][5])*cd14*flags->swc[5];
        t82 = (p[33]*ctx->plg[2][3]+p[36]*ctx->plg[2][5])*cd14*flags->swc[5];
        t[7] = ((p[5]*ctx->plg[2][2] + p[41]*ctx->plg[2][4] + t81) * ctx->c2tloc + (p[8]*ctx->plg[2][2] + p[42]*ctx->plg[2][4] + t82) * ctx->s2tloc);
    }

    /* TERDIURNAL */
    if (flags->sw[14]) {
        t[13] = p[39] * ctx->plg[3][3] * ctx->s3tloc + p[40] * ctx->plg[3][3] * ctx->c3tloc;
    }

    /* MAGNETIC ACTIVITY */
    if (flags->sw[9]) {
        if (flags->sw[9]==1)
            t[8] = ctx->apdf * (p[32] + p[45] * ctx->plg[0][2] * flags->swc[2]);
        if (flags->sw[9]==-1)
            t[8]=(p[50]*ctx->apt[0] + p[96]*ctx->plg[0][2] * ctx->apt[0]*flags->swc[2]);
    }

    /* LONGITUDINAL */
    if (!((flags->sw[10]==0) || (flags->sw[11]==0) || (input->g_long<=-1000.0))) {
        t[10] = (1.0 + ctx->plg[0][1]*(p[80]*flags->swc[5]*cos(dr*(input->doy-p[81]))\
                +p[85]*flags->swc[6]*cos(2.0*dr*(input->doy-p[86])))\
            +p[83]*flags->swc[3]*cos(dr*(input->doy-p[84]))\
            +p[87]*flags->swc[4]*cos(2.0*dr*(input->doy-p[88])))\
            *((p[64]*ctx->plg[1][2]+p[65]*ctx->plg[1][4]+p[66]*ctx->plg[1][6]\
            +p[74]*ctx->plg[1][1]+p[75]*ctx->plg[1][3]+p[76]*ctx->plg[1][5]\
            )*cos(dgtr*input->g_long)\
            +(p[90]*ctx->plg[1][2]+p[91]*ctx->plg[1][4]+p[92]*ctx->plg[1][6]\
            +p[77]*ctx->plg[1][1]+p[78]*ctx->plg[1][3]+p[79]*ctx->plg[1][5]\
            )*sin(dgtr*input->g_long));
    }
    tt=0;
    for (i=0;i<14;i++)
        tt+=fabs(flags->sw[i+1])*t[i];
    return tt;
}



/* ------------------------------------------------------------------- */
/* ------------------------------- GTD7 ------------------------------ */
/* ------------------------------------------------------------------- */

void tinymsis_core_gtd7(tinymsis_workspace *ctx, struct nrlmsise_input *input, struct nrlmsise_flags *flags, struct nrlmsise_output *output) {
    (void)ctx;
    double xlat;
    double xmm;
    int mn3 = 5;
    double zn3[5]={32.5,20.0,15.0,10.0,0.0};
    int mn2 = 4;
    double zn2[4]={72.5,55.0,45.0,32.5};
    double altt;
    double zmix=62.5;
    double tmp;
    double dm28m;
    double tz;
    double dmc;
    double dmr;
    double dz28;
    struct nrlmsise_output soutput;
    int i;

    tinymsis_core_tselec(ctx, flags);

    /* Latitude variation of gravity (none for sw[2]=0) */
    xlat=input->g_lat;
    if (flags->sw[2]==0)
        xlat=45.0;
    tinymsis_core_glatf(ctx, xlat, &ctx->gsurf, &ctx->re);

    xmm = tinymsis_coeff_pdm[2][4];

    /* THERMOSPHERE / MESOSPHERE (above zn2[0]) */
    if (input->alt>zn2[0])
        altt=input->alt;
    else
        altt=zn2[0];

    tmp=input->alt;
    input->alt=altt;
    tinymsis_core_gts7(ctx, input, flags, &soutput);
    altt=input->alt;
    input->alt=tmp;
    if (flags->sw[0])   /* metric adjustment */
        dm28m=ctx->dm28*1.0E6;
    else
        dm28m=ctx->dm28;
    output->t[0]=soutput.t[0];
    output->t[1]=soutput.t[1];
    if (input->alt>=zn2[0]) {
        for (i=0;i<9;i++)
            output->d[i]=soutput.d[i];
        return;
    }

/*       LOWER MESOSPHERE/UPPER STRATOSPHERE (between zn3[0] and zn2[0])
 *         Temperature at nodes and gradients at end nodes
 *         Inverse temperature a linear function of spherical harmonics
 */
    ctx->meso_tgn2[0]=ctx->meso_tgn1[1];
    ctx->meso_tn2[0]=ctx->meso_tn1[4];
        ctx->meso_tn2[1]=tinymsis_coeff_pma[0][0]*tinymsis_coeff_pavgm[0]/(1.0-flags->sw[20]*tinymsis_core_glob7s(ctx, tinymsis_coeff_pma[0], input, flags));
        ctx->meso_tn2[2]=tinymsis_coeff_pma[1][0]*tinymsis_coeff_pavgm[1]/(1.0-flags->sw[20]*tinymsis_core_glob7s(ctx, tinymsis_coeff_pma[1], input, flags));
        ctx->meso_tn2[3]=tinymsis_coeff_pma[2][0]*tinymsis_coeff_pavgm[2]/(1.0-flags->sw[20]*flags->sw[22]*tinymsis_core_glob7s(ctx, tinymsis_coeff_pma[2], input, flags));
    ctx->meso_tgn2[1]=tinymsis_coeff_pavgm[8]*tinymsis_coeff_pma[9][0]*(1.0+flags->sw[20]*flags->sw[22]*tinymsis_core_glob7s(ctx, tinymsis_coeff_pma[9], input, flags))*ctx->meso_tn2[3]*ctx->meso_tn2[3]/(pow((tinymsis_coeff_pma[2][0]*tinymsis_coeff_pavgm[2]),2.0));
    ctx->meso_tn3[0]=ctx->meso_tn2[3];

    if (input->alt<=zn3[0]) {
/*       LOWER STRATOSPHERE AND TROPOSPHERE (below zn3[0])
 *         Temperature at nodes and gradients at end nodes
 *         Inverse temperature a linear function of spherical harmonics
 */
        ctx->meso_tgn3[0]=ctx->meso_tgn2[1];
        ctx->meso_tn3[1]=tinymsis_coeff_pma[3][0]*tinymsis_coeff_pavgm[3]/(1.0-flags->sw[22]*tinymsis_core_glob7s(ctx, tinymsis_coeff_pma[3], input, flags));
        ctx->meso_tn3[2]=tinymsis_coeff_pma[4][0]*tinymsis_coeff_pavgm[4]/(1.0-flags->sw[22]*tinymsis_core_glob7s(ctx, tinymsis_coeff_pma[4], input, flags));
        ctx->meso_tn3[3]=tinymsis_coeff_pma[5][0]*tinymsis_coeff_pavgm[5]/(1.0-flags->sw[22]*tinymsis_core_glob7s(ctx, tinymsis_coeff_pma[5], input, flags));
        ctx->meso_tn3[4]=tinymsis_coeff_pma[6][0]*tinymsis_coeff_pavgm[6]/(1.0-flags->sw[22]*tinymsis_core_glob7s(ctx, tinymsis_coeff_pma[6], input, flags));
        ctx->meso_tgn3[1]=tinymsis_coeff_pma[7][0]*tinymsis_coeff_pavgm[7]*(1.0+flags->sw[22]*tinymsis_core_glob7s(ctx, tinymsis_coeff_pma[7], input, flags)) *ctx->meso_tn3[4]*ctx->meso_tn3[4]/(pow((tinymsis_coeff_pma[6][0]*tinymsis_coeff_pavgm[6]),2.0));
    }

        /* LINEAR TRANSITION TO FULL MIXING BELOW zn2[0] */

    dmc=0;
    if (input->alt>zmix)
        dmc = 1.0 - (zn2[0]-input->alt)/(zn2[0] - zmix);
    dz28=soutput.d[2];

    /**** N2 density ****/
    dmr=soutput.d[2] / dm28m - 1.0;
    output->d[2]=tinymsis_core_densm(ctx, input->alt,dm28m,xmm, &tz, mn3, zn3, ctx->meso_tn3, ctx->meso_tgn3, mn2, zn2, ctx->meso_tn2, ctx->meso_tgn2);
    output->d[2]=output->d[2] * (1.0 + dmr*dmc);

    /**** HE density ****/
    dmr = soutput.d[0] / (dz28 * tinymsis_coeff_pdm[0][1]) - 1.0;
    output->d[0] = output->d[2] * tinymsis_coeff_pdm[0][1] * (1.0 + dmr*dmc);

    /**** O density ****/
    output->d[1] = 0;
    output->d[8] = 0;

    /**** O2 density ****/
    dmr = soutput.d[3] / (dz28 * tinymsis_coeff_pdm[3][1]) - 1.0;
    output->d[3] = output->d[2] * tinymsis_coeff_pdm[3][1] * (1.0 + dmr*dmc);

    /**** AR density ***/
    dmr = soutput.d[4] / (dz28 * tinymsis_coeff_pdm[4][1]) - 1.0;
    output->d[4] = output->d[2] * tinymsis_coeff_pdm[4][1] * (1.0 + dmr*dmc);

    /**** Hydrogen density ****/
    output->d[6] = 0;

    /**** Atomic nitrogen density ****/
    output->d[7] = 0;

    /**** Total mass density */
    output->d[5] = 1.66E-24 * (4.0 * output->d[0] + 16.0 * output->d[1] + 28.0 * output->d[2] + 32.0 * output->d[3] + 40.0 * output->d[4] + output->d[6] + 14.0 * output->d[7]);

    if (flags->sw[0])
        output->d[5]=output->d[5]/1000;

    /**** temperature at altitude ****/
    ctx->dd = tinymsis_core_densm(ctx, input->alt, 1.0, 0, &tz, mn3, zn3, ctx->meso_tn3, ctx->meso_tgn3, mn2, zn2, ctx->meso_tn2, ctx->meso_tgn2);
    output->t[1]=tz;

}



/* ------------------------------------------------------------------- */
/* ------------------------------- GTD7D ----------------------------- */
/* ------------------------------------------------------------------- */

void tinymsis_core_gtd7d(tinymsis_workspace *ctx, struct nrlmsise_input *input, struct nrlmsise_flags *flags, struct nrlmsise_output *output) {
    (void)ctx;
    tinymsis_core_gtd7(ctx, input, flags, output);
    output->d[5] = 1.66E-24 * (4.0 * output->d[0] + 16.0 * output->d[1] + 28.0 * output->d[2] + 32.0 * output->d[3] + 40.0 * output->d[4] + output->d[6] + 14.0 * output->d[7] + 16.0 * output->d[8]);
    if (flags->sw[0])
        output->d[5]=output->d[5]/1000;
}



/* ------------------------------------------------------------------- */
/* -------------------------------- GHP7 ----------------------------- */
/* ------------------------------------------------------------------- */





/* ------------------------------------------------------------------- */
/* ------------------------------- GTS7 ------------------------------ */
/* ------------------------------------------------------------------- */

static void tinymsis_core_gts7(tinymsis_workspace *ctx, struct nrlmsise_input *input, struct nrlmsise_flags *flags, struct nrlmsise_output *output) {
    (void)ctx;
/*     Thermospheric portion of NRLMSISE-00
 *     See GTD7 for more extensive comments
 *     alt > 72.5 km!
 */
    double za;
    int i, j;
    double ddum, z;
    double zn1[5] = {120.0, 110.0, 100.0, 90.0, 72.5};
    double tinf;
    int mn1 = 5;
    double tinymsis_core_g0;
    double tlb;
    double s;
    double db01, db04, db14, db16, db28, db32, db40;
    double zh28, zh04, zh16, zh32, zh40, zh01, zh14;
    double zhm28, zhm04, zhm16, zhm32, zhm40, zhm01, zhm14;
    double xmd;
    double b28, b04, b16, b32, b40, b01, b14;
    double tz;
    double g28, g4, g16, g32, g40, g1, g14;
    double zhf, xmm;
    double zc04, zc16, zc32, zc40, zc01, zc14;
    double hc04, hc16, hc32, hc40, hc01, hc14;
    double hcc16, hcc32, hcc01, hcc14;
    double zcc16, zcc32, zcc01, zcc14;
    double rc16, rc32, rc01, rc14;
    double rl;
    double g16h, db16h, tho, zsht, zmho, zsho;
    double dgtr=1.74533E-2;
    double dr=1.72142E-2;
    double alpha[9]={-0.38, 0.0, 0.0, 0.0, 0.17, 0.0, -0.38, 0.0, 0.0};
    double altl[8]={200.0, 300.0, 160.0, 250.0, 240.0, 450.0, 320.0, 450.0};
    double diffusive_density;
    double hc216, hcc232;
    za = tinymsis_coeff_pdl[1][15];
    zn1[0] = za;
    for (j=0;j<9;j++)
        output->d[j]=0;

    /* TINF VARIATIONS NOT IMPORTANT BELOW ZA OR ZN1(1) */
    if (input->alt>zn1[0])
        tinf = tinymsis_coeff_ptm[0]*tinymsis_coeff_pt[0] * \
            (1.0+flags->sw[16]*tinymsis_core_globe7(ctx, tinymsis_coeff_pt,input,flags));
    else
        tinf = tinymsis_coeff_ptm[0]*tinymsis_coeff_pt[0];
    output->t[0]=tinf;

    /*  GRADIENT VARIATIONS NOT IMPORTANT BELOW ZN1(5) */
    if (input->alt>zn1[4])
        tinymsis_core_g0 = tinymsis_coeff_ptm[3]*tinymsis_coeff_ps[0] * \
            (1.0+flags->sw[19]*tinymsis_core_globe7(ctx, tinymsis_coeff_ps,input,flags));
    else
        tinymsis_core_g0 = tinymsis_coeff_ptm[3]*tinymsis_coeff_ps[0];
    tlb = tinymsis_coeff_ptm[1] * (1.0 + flags->sw[17]*tinymsis_core_globe7(ctx, tinymsis_coeff_pd[3],input,flags))*tinymsis_coeff_pd[3][0];
    s = tinymsis_core_g0 / (tinf - tlb);

/*      Lower thermosphere temp variations not significant for
 *       density above 300 km */
    if (input->alt<300.0) {
        ctx->meso_tn1[1]=tinymsis_coeff_ptm[6]*tinymsis_coeff_ptl[0][0]/(1.0-flags->sw[18]*tinymsis_core_glob7s(ctx, tinymsis_coeff_ptl[0], input, flags));
        ctx->meso_tn1[2]=tinymsis_coeff_ptm[2]*tinymsis_coeff_ptl[1][0]/(1.0-flags->sw[18]*tinymsis_core_glob7s(ctx, tinymsis_coeff_ptl[1], input, flags));
        ctx->meso_tn1[3]=tinymsis_coeff_ptm[7]*tinymsis_coeff_ptl[2][0]/(1.0-flags->sw[18]*tinymsis_core_glob7s(ctx, tinymsis_coeff_ptl[2], input, flags));
        ctx->meso_tn1[4]=tinymsis_coeff_ptm[4]*tinymsis_coeff_ptl[3][0]/(1.0-flags->sw[18]*flags->sw[20]*tinymsis_core_glob7s(ctx, tinymsis_coeff_ptl[3], input, flags));
        ctx->meso_tgn1[1]=tinymsis_coeff_ptm[8]*tinymsis_coeff_pma[8][0]*(1.0+flags->sw[18]*flags->sw[20]*tinymsis_core_glob7s(ctx, tinymsis_coeff_pma[8], input, flags))*ctx->meso_tn1[4]*ctx->meso_tn1[4]/(pow((tinymsis_coeff_ptm[4]*tinymsis_coeff_ptl[3][0]),2.0));
    } else {
        ctx->meso_tn1[1]=tinymsis_coeff_ptm[6]*tinymsis_coeff_ptl[0][0];
        ctx->meso_tn1[2]=tinymsis_coeff_ptm[2]*tinymsis_coeff_ptl[1][0];
        ctx->meso_tn1[3]=tinymsis_coeff_ptm[7]*tinymsis_coeff_ptl[2][0];
        ctx->meso_tn1[4]=tinymsis_coeff_ptm[4]*tinymsis_coeff_ptl[3][0];
        ctx->meso_tgn1[1]=tinymsis_coeff_ptm[8]*tinymsis_coeff_pma[8][0]*ctx->meso_tn1[4]*ctx->meso_tn1[4]/(pow((tinymsis_coeff_ptm[4]*tinymsis_coeff_ptl[3][0]),2.0));
    }

    /* N2 variation factor at Zlb */
    g28=flags->sw[21]*tinymsis_core_globe7(ctx, tinymsis_coeff_pd[2], input, flags);

    /* VARIATION OF TURBOPAUSE HEIGHT */
    zhf=tinymsis_coeff_pdl[1][24]*(1.0+flags->sw[5]*tinymsis_coeff_pdl[0][24]*sin(dgtr*input->g_lat)*cos(dr*(input->doy-tinymsis_coeff_pt[13])));
    output->t[0]=tinf;
    xmm = tinymsis_coeff_pdm[2][4];
    z = input->alt;


        /**** N2 DENSITY ****/

    /* Diffusive density at Zlb */
    db28 = tinymsis_coeff_pdm[2][0]*exp(g28)*tinymsis_coeff_pd[2][0];
    /* Diffusive density at Alt */
    output->d[2]=tinymsis_core_densu(ctx, z,db28,tinf,tlb,28.0,alpha[2],&output->t[1],tinymsis_coeff_ptm[5],s,mn1,zn1,ctx->meso_tn1,ctx->meso_tgn1);
    diffusive_density=output->d[2];
    /* Turbopause */
    zh28=tinymsis_coeff_pdm[2][2]*zhf;
    zhm28=tinymsis_coeff_pdm[2][3]*tinymsis_coeff_pdl[1][5];
    xmd=28.0-xmm;
    /* Mixed density at Zlb */
    b28=tinymsis_core_densu(ctx, zh28,db28,tinf,tlb,xmd,(alpha[2]-1.0),&tz,tinymsis_coeff_ptm[5],s,mn1, zn1,ctx->meso_tn1,ctx->meso_tgn1);
    if ((flags->sw[15])&&(z<=altl[2])) {
        /*  Mixed density at Alt */
        ctx->dm28=tinymsis_core_densu(ctx, z,b28,tinf,tlb,xmm,alpha[2],&tz,tinymsis_coeff_ptm[5],s,mn1,zn1,ctx->meso_tn1,ctx->meso_tgn1);
        /*  Net density at Alt */
        output->d[2]=tinymsis_core_dnet(ctx, output->d[2],ctx->dm28,zhm28,xmm,28.0);
    }


        /**** HE DENSITY ****/

    /*   Density variation factor at Zlb */
    g4 = flags->sw[21]*tinymsis_core_globe7(ctx, tinymsis_coeff_pd[0], input, flags);
    /*  Diffusive density at Zlb */
    db04 = tinymsis_coeff_pdm[0][0]*exp(g4)*tinymsis_coeff_pd[0][0];
        /*  Diffusive density at Alt */
    output->d[0]=tinymsis_core_densu(ctx, z,db04,tinf,tlb, 4.,alpha[0],&output->t[1],tinymsis_coeff_ptm[5],s,mn1,zn1,ctx->meso_tn1,ctx->meso_tgn1);
    diffusive_density=output->d[0];
    if ((flags->sw[15]) && (z<altl[0])) {
        /*  Turbopause */
        zh04=tinymsis_coeff_pdm[0][2];
        /*  Mixed density at Zlb */
        b04=tinymsis_core_densu(ctx, zh04,db04,tinf,tlb,4.-xmm,alpha[0]-1.,&output->t[1],tinymsis_coeff_ptm[5],s,mn1,zn1,ctx->meso_tn1,ctx->meso_tgn1);
        /*  Mixed density at Alt */
        ctx->dm04=tinymsis_core_densu(ctx, z,b04,tinf,tlb,xmm,0.,&output->t[1],tinymsis_coeff_ptm[5],s,mn1,zn1,ctx->meso_tn1,ctx->meso_tgn1);
        zhm04=zhm28;
        /*  Net density at Alt */
        output->d[0]=tinymsis_core_dnet(ctx, output->d[0],ctx->dm04,zhm04,xmm,4.);
        /*  Correction to specified mixing ratio at ground */
        rl=log(b28*tinymsis_coeff_pdm[0][1]/b04);
        zc04=tinymsis_coeff_pdm[0][4]*tinymsis_coeff_pdl[1][0];
        hc04=tinymsis_coeff_pdm[0][5]*tinymsis_coeff_pdl[1][1];
        /*  Net density corrected at Alt */
        output->d[0]=output->d[0]*tinymsis_core_ccor(ctx, z,rl,hc04,zc04);
    }


        /**** O DENSITY ****/

    /*  Density variation factor at Zlb */
    g16= flags->sw[21]*tinymsis_core_globe7(ctx, tinymsis_coeff_pd[1],input,flags);
    /*  Diffusive density at Zlb */
    db16 =  tinymsis_coeff_pdm[1][0]*exp(g16)*tinymsis_coeff_pd[1][0];
        /*   Diffusive density at Alt */
    output->d[1]=tinymsis_core_densu(ctx, z,db16,tinf,tlb, 16.,alpha[1],&output->t[1],tinymsis_coeff_ptm[5],s,mn1, zn1,ctx->meso_tn1,ctx->meso_tgn1);
    diffusive_density=output->d[1];
    if ((flags->sw[15]) && (z<=altl[1])) {
        /*   Turbopause */
        zh16=tinymsis_coeff_pdm[1][2];
        /*  Mixed density at Zlb */
        b16=tinymsis_core_densu(ctx, zh16,db16,tinf,tlb,16.0-xmm,(alpha[1]-1.0), &output->t[1],tinymsis_coeff_ptm[5],s,mn1,zn1,ctx->meso_tn1,ctx->meso_tgn1);
        /*  Mixed density at Alt */
        ctx->dm16=tinymsis_core_densu(ctx, z,b16,tinf,tlb,xmm,0.,&output->t[1],tinymsis_coeff_ptm[5],s,mn1,zn1,ctx->meso_tn1,ctx->meso_tgn1);
        zhm16=zhm28;
        /*  Net density at Alt */
        output->d[1]=tinymsis_core_dnet(ctx, output->d[1],ctx->dm16,zhm16,xmm,16.);
        rl=tinymsis_coeff_pdm[1][1]*tinymsis_coeff_pdl[1][16]*(1.0+flags->sw[1]*tinymsis_coeff_pdl[0][23]*(input->f107A-150.0));
        hc16=tinymsis_coeff_pdm[1][5]*tinymsis_coeff_pdl[1][3];
        zc16=tinymsis_coeff_pdm[1][4]*tinymsis_coeff_pdl[1][2];
        hc216=tinymsis_coeff_pdm[1][5]*tinymsis_coeff_pdl[1][4];
        output->d[1]=output->d[1]*tinymsis_core_ccor2(ctx, z,rl,hc16,zc16,hc216);
        /*   Chemistry correction */
        hcc16=tinymsis_coeff_pdm[1][7]*tinymsis_coeff_pdl[1][13];
        zcc16=tinymsis_coeff_pdm[1][6]*tinymsis_coeff_pdl[1][12];
        rc16=tinymsis_coeff_pdm[1][3]*tinymsis_coeff_pdl[1][14];
        /*  Net density corrected at Alt */
        output->d[1]=output->d[1]*tinymsis_core_ccor(ctx, z,rc16,hcc16,zcc16);
    }


        /**** O2 DENSITY ****/

        /*   Density variation factor at Zlb */
    g32= flags->sw[21]*tinymsis_core_globe7(ctx, tinymsis_coeff_pd[4], input, flags);
        /*  Diffusive density at Zlb */
    db32 = tinymsis_coeff_pdm[3][0]*exp(g32)*tinymsis_coeff_pd[4][0];
        /*   Diffusive density at Alt */
    output->d[3]=tinymsis_core_densu(ctx, z,db32,tinf,tlb, 32.,alpha[3],&output->t[1],tinymsis_coeff_ptm[5],s,mn1, zn1,ctx->meso_tn1,ctx->meso_tgn1);
    diffusive_density=output->d[3];
    if (flags->sw[15]) {
        if (z<=altl[3]) {
            /*   Turbopause */
            zh32=tinymsis_coeff_pdm[3][2];
            /*  Mixed density at Zlb */
            b32=tinymsis_core_densu(ctx, zh32,db32,tinf,tlb,32.-xmm,alpha[3]-1., &output->t[1],tinymsis_coeff_ptm[5],s,mn1,zn1,ctx->meso_tn1,ctx->meso_tgn1);
            /*  Mixed density at Alt */
            ctx->dm32=tinymsis_core_densu(ctx, z,b32,tinf,tlb,xmm,0.,&output->t[1],tinymsis_coeff_ptm[5],s,mn1,zn1,ctx->meso_tn1,ctx->meso_tgn1);
            zhm32=zhm28;
            /*  Net density at Alt */
            output->d[3]=tinymsis_core_dnet(ctx, output->d[3],ctx->dm32,zhm32,xmm,32.);
            /*   Correction to specified mixing ratio at ground */
            rl=log(b28*tinymsis_coeff_pdm[3][1]/b32);
            hc32=tinymsis_coeff_pdm[3][5]*tinymsis_coeff_pdl[1][7];
            zc32=tinymsis_coeff_pdm[3][4]*tinymsis_coeff_pdl[1][6];
            output->d[3]=output->d[3]*tinymsis_core_ccor(ctx, z,rl,hc32,zc32);
        }
        /*  Correction for general departure from diffusive equilibrium above Zlb */
        hcc32=tinymsis_coeff_pdm[3][7]*tinymsis_coeff_pdl[1][22];
        hcc232=tinymsis_coeff_pdm[3][7]*tinymsis_coeff_pdl[0][22];
        zcc32=tinymsis_coeff_pdm[3][6]*tinymsis_coeff_pdl[1][21];
        rc32=tinymsis_coeff_pdm[3][3]*tinymsis_coeff_pdl[1][23]*(1.+flags->sw[1]*tinymsis_coeff_pdl[0][23]*(input->f107A-150.));
        /*  Net density corrected at Alt */
        output->d[3]=output->d[3]*tinymsis_core_ccor2(ctx, z,rc32,hcc32,zcc32,hcc232);
    }


        /**** AR DENSITY ****/

        /*   Density variation factor at Zlb */
    g40= flags->sw[21]*tinymsis_core_globe7(ctx, tinymsis_coeff_pd[5],input,flags);
        /*  Diffusive density at Zlb */
    db40 = tinymsis_coeff_pdm[4][0]*exp(g40)*tinymsis_coeff_pd[5][0];
    /*   Diffusive density at Alt */
    output->d[4]=tinymsis_core_densu(ctx, z,db40,tinf,tlb, 40.,alpha[4],&output->t[1],tinymsis_coeff_ptm[5],s,mn1,zn1,ctx->meso_tn1,ctx->meso_tgn1);
    diffusive_density=output->d[4];
    if ((flags->sw[15]) && (z<=altl[4])) {
        /*   Turbopause */
        zh40=tinymsis_coeff_pdm[4][2];
        /*  Mixed density at Zlb */
        b40=tinymsis_core_densu(ctx, zh40,db40,tinf,tlb,40.-xmm,alpha[4]-1.,&output->t[1],tinymsis_coeff_ptm[5],s,mn1,zn1,ctx->meso_tn1,ctx->meso_tgn1);
        /*  Mixed density at Alt */
        ctx->dm40=tinymsis_core_densu(ctx, z,b40,tinf,tlb,xmm,0.,&output->t[1],tinymsis_coeff_ptm[5],s,mn1,zn1,ctx->meso_tn1,ctx->meso_tgn1);
        zhm40=zhm28;
        /*  Net density at Alt */
        output->d[4]=tinymsis_core_dnet(ctx, output->d[4],ctx->dm40,zhm40,xmm,40.);
        /*   Correction to specified mixing ratio at ground */
        rl=log(b28*tinymsis_coeff_pdm[4][1]/b40);
        hc40=tinymsis_coeff_pdm[4][5]*tinymsis_coeff_pdl[1][9];
        zc40=tinymsis_coeff_pdm[4][4]*tinymsis_coeff_pdl[1][8];
        /*  Net density corrected at Alt */
        output->d[4]=output->d[4]*tinymsis_core_ccor(ctx, z,rl,hc40,zc40);
      }


        /**** HYDROGEN DENSITY ****/

        /*   Density variation factor at Zlb */
    g1 = flags->sw[21]*tinymsis_core_globe7(ctx, tinymsis_coeff_pd[6], input, flags);
        /*  Diffusive density at Zlb */
    db01 = tinymsis_coeff_pdm[5][0]*exp(g1)*tinymsis_coeff_pd[6][0];
        /*   Diffusive density at Alt */
    output->d[6]=tinymsis_core_densu(ctx, z,db01,tinf,tlb,1.,alpha[6],&output->t[1],tinymsis_coeff_ptm[5],s,mn1,zn1,ctx->meso_tn1,ctx->meso_tgn1);
    diffusive_density=output->d[6];
    if ((flags->sw[15]) && (z<=altl[6])) {
        /*   Turbopause */
        zh01=tinymsis_coeff_pdm[5][2];
        /*  Mixed density at Zlb */
        b01=tinymsis_core_densu(ctx, zh01,db01,tinf,tlb,1.-xmm,alpha[6]-1., &output->t[1],tinymsis_coeff_ptm[5],s,mn1,zn1,ctx->meso_tn1,ctx->meso_tgn1);
        /*  Mixed density at Alt */
        ctx->dm01=tinymsis_core_densu(ctx, z,b01,tinf,tlb,xmm,0.,&output->t[1],tinymsis_coeff_ptm[5],s,mn1,zn1,ctx->meso_tn1,ctx->meso_tgn1);
        zhm01=zhm28;
        /*  Net density at Alt */
        output->d[6]=tinymsis_core_dnet(ctx, output->d[6],ctx->dm01,zhm01,xmm,1.);
        /*   Correction to specified mixing ratio at ground */
        rl=log(b28*tinymsis_coeff_pdm[5][1]*sqrt(tinymsis_coeff_pdl[1][17]*tinymsis_coeff_pdl[1][17])/b01);
        hc01=tinymsis_coeff_pdm[5][5]*tinymsis_coeff_pdl[1][11];
        zc01=tinymsis_coeff_pdm[5][4]*tinymsis_coeff_pdl[1][10];
        output->d[6]=output->d[6]*tinymsis_core_ccor(ctx, z,rl,hc01,zc01);
        /*   Chemistry correction */
        hcc01=tinymsis_coeff_pdm[5][7]*tinymsis_coeff_pdl[1][19];
        zcc01=tinymsis_coeff_pdm[5][6]*tinymsis_coeff_pdl[1][18];
        rc01=tinymsis_coeff_pdm[5][3]*tinymsis_coeff_pdl[1][20];
        /*  Net density corrected at Alt */
        output->d[6]=output->d[6]*tinymsis_core_ccor(ctx, z,rc01,hcc01,zcc01);
}


        /**** ATOMIC NITROGEN DENSITY ****/

    /*   Density variation factor at Zlb */
    g14 = flags->sw[21]*tinymsis_core_globe7(ctx, tinymsis_coeff_pd[7],input,flags);
        /*  Diffusive density at Zlb */
    db14 = tinymsis_coeff_pdm[6][0]*exp(g14)*tinymsis_coeff_pd[7][0];
        /*   Diffusive density at Alt */
    output->d[7]=tinymsis_core_densu(ctx, z,db14,tinf,tlb,14.,alpha[7],&output->t[1],tinymsis_coeff_ptm[5],s,mn1,zn1,ctx->meso_tn1,ctx->meso_tgn1);
    diffusive_density=output->d[7];
    if ((flags->sw[15]) && (z<=altl[7])) {
        /*   Turbopause */
        zh14=tinymsis_coeff_pdm[6][2];
        /*  Mixed density at Zlb */
        b14=tinymsis_core_densu(ctx, zh14,db14,tinf,tlb,14.-xmm,alpha[7]-1., &output->t[1],tinymsis_coeff_ptm[5],s,mn1,zn1,ctx->meso_tn1,ctx->meso_tgn1);
        /*  Mixed density at Alt */
        ctx->dm14=tinymsis_core_densu(ctx, z,b14,tinf,tlb,xmm,0.,&output->t[1],tinymsis_coeff_ptm[5],s,mn1,zn1,ctx->meso_tn1,ctx->meso_tgn1);
        zhm14=zhm28;
        /*  Net density at Alt */
        output->d[7]=tinymsis_core_dnet(ctx, output->d[7],ctx->dm14,zhm14,xmm,14.);
        /*   Correction to specified mixing ratio at ground */
        rl=log(b28*tinymsis_coeff_pdm[6][1]*sqrt(tinymsis_coeff_pdl[0][2]*tinymsis_coeff_pdl[0][2])/b14);
        hc14=tinymsis_coeff_pdm[6][5]*tinymsis_coeff_pdl[0][1];
        zc14=tinymsis_coeff_pdm[6][4]*tinymsis_coeff_pdl[0][0];
        output->d[7]=output->d[7]*tinymsis_core_ccor(ctx, z,rl,hc14,zc14);
        /*   Chemistry correction */
        hcc14=tinymsis_coeff_pdm[6][7]*tinymsis_coeff_pdl[0][4];
        zcc14=tinymsis_coeff_pdm[6][6]*tinymsis_coeff_pdl[0][3];
        rc14=tinymsis_coeff_pdm[6][3]*tinymsis_coeff_pdl[0][5];
        /*  Net density corrected at Alt */
        output->d[7]=output->d[7]*tinymsis_core_ccor(ctx, z,rc14,hcc14,zcc14);
    }


        /**** Anomalous OXYGEN DENSITY ****/

    g16h = flags->sw[21]*tinymsis_core_globe7(ctx, tinymsis_coeff_pd[8],input,flags);
    db16h = tinymsis_coeff_pdm[7][0]*exp(g16h)*tinymsis_coeff_pd[8][0];
    tho = tinymsis_coeff_pdm[7][9]*tinymsis_coeff_pdl[0][6];
    diffusive_density=tinymsis_core_densu(ctx, z,db16h,tho,tho,16.,alpha[8],&output->t[1],tinymsis_coeff_ptm[5],s,mn1, zn1,ctx->meso_tn1,ctx->meso_tgn1);
    zsht=tinymsis_coeff_pdm[7][5];
    zmho=tinymsis_coeff_pdm[7][4];
    zsho=tinymsis_core_scalh(ctx, zmho,16.0,tho);
    output->d[8]=diffusive_density*exp(-zsht/zsho*(exp(-(z-zmho)/zsht)-1.));


    /* total mass density */
    output->d[5] = 1.66E-24*(4.0*output->d[0]+16.0*output->d[1]+28.0*output->d[2]+32.0*output->d[3]+40.0*output->d[4]+ output->d[6]+14.0*output->d[7]);


    /* temperature */
    z = sqrt(input->alt*input->alt);
    ddum = tinymsis_core_densu(ctx, z,1.0, tinf, tlb, 0.0, 0.0, &output->t[1], tinymsis_coeff_ptm[5], s, mn1, zn1, ctx->meso_tn1, ctx->meso_tgn1);
    (void) ddum; /* silence gcc */
    if (flags->sw[0]) {
        for(i=0;i<9;i++)
            output->d[i]=output->d[i]*1.0E6;
        output->d[5]=output->d[5]/1000;
    }
}
