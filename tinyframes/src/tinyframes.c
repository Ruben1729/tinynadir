/* SPDX-License-Identifier: MIT; private ERFA subset has separate notices. */
#include "tinyframes/tinyframes.h"
#include "erfa.h"
#include <math.h>
#include <stddef.h>
static void normalize(double in0,double in1,double out[2]) {
    double a=floor(in0-.5)+.5,b=(in0-a)+in1,w=floor(b);out[0]=a+w;out[1]=b-w;
}
static tinyframes_status valid(const tinyframes_epoch *in,tinyframes_frame from,tinyframes_frame to,tinyframes_model model,tinyframes_epoch *e) {
    double t,u,d;int i;
    if(!in)return TINYFRAMES_INVALID;
    if((unsigned)from>(unsigned)TINYFRAMES_PEF1982||(unsigned)to>(unsigned)TINYFRAMES_PEF1982||(model!=TINYFRAMES_IAU2006A&&model!=TINYFRAMES_IAU2000B))return TINYFRAMES_RANGE;
    for(i=0;i<2;i++)if(!isfinite(in->tt[i])||!isfinite(in->ut1[i]))return TINYFRAMES_INVALID;
    if(!isfinite(in->xp_rad)||!isfinite(in->yp_rad)||!isfinite(in->dx_rad)||!isfinite(in->dy_rad)||!isfinite(in->lod_s))return TINYFRAMES_INVALID;
    if(fabs(in->tt[0])>1e8||fabs(in->tt[1])>1e8||fabs(in->ut1[0])>1e8||fabs(in->ut1[1])>1e8)return TINYFRAMES_RANGE;
    t=((in->tt[0]-2451545)+in->tt[1])/36525;u=((in->ut1[0]-2451545)+in->ut1[1])/36525;
    d=(in->tt[0]-in->ut1[0])+(in->tt[1]-in->ut1[1]);
    if(fabs(t)>1||fabs(u)>1||fabs(d)>1||fabs(in->xp_rad)>.01||fabs(in->yp_rad)>.01||fabs(in->dx_rad)>.001||fabs(in->dy_rad)>.001||fabs(in->lod_s)>10)return TINYFRAMES_RANGE;
    *e=*in;normalize(in->tt[0],in->tt[1],e->tt);normalize(in->ut1[0],in->ut1[1],e->ut1);return TINYFRAMES_OK;
}
static void copy(double a[3][3],double out[3][3]) {int i,j;for(i=0;i<3;i++)for(j=0;j<3;j++)out[i][j]=a[i][j];}
static void multiply(double a[3][3],double b[3][3],double out[3][3]) {
    double r[3][3];int i,j,k;for(i=0;i<3;i++)for(j=0;j<3;j++){r[i][j]=0;for(k=0;k<3;k++)r[i][j]+=a[i][k]*b[k][j];}copy(r,out);
}
static void transpose(double a[3][3],double out[3][3]) {int i,j;for(i=0;i<3;i++)for(j=0;j<3;j++)out[i][j]=a[j][i];}
static void model_matrices(const tinyframes_epoch *e,tinyframes_model model,double ci[3][3],double pom[3][3],double gcrf[3][3]) {
    double bpn[3][3],x,y,s,era;
    if(model==TINYFRAMES_IAU2006A)tinyframes_erfa_Pnm06a(e->tt[0],e->tt[1],bpn);
    else tinyframes_erfa_Pnm00b(e->tt[0],e->tt[1],bpn);
    tinyframes_erfa_Bpn2xy(bpn,&x,&y);x+=e->dx_rad;y+=e->dy_rad;
    s=model==TINYFRAMES_IAU2006A?tinyframes_erfa_S06(e->tt[0],e->tt[1],x,y):tinyframes_erfa_S00(e->tt[0],e->tt[1],x,y);
    tinyframes_erfa_C2ixys(x,y,s,ci);
    tinyframes_erfa_Pom00(e->xp_rad,e->yp_rad,model==TINYFRAMES_IAU2006A?tinyframes_erfa_Sp00(e->tt[0],e->tt[1]):0,pom);
    era=tinyframes_erfa_Era00(e->ut1[0],e->ut1[1]);
    tinyframes_erfa_C2tcio(ci,era,pom,gcrf);
}
static void to_fixed(const tinyframes_epoch *e,tinyframes_frame frame,tinyframes_model model,double out[3][3]) {
    double ci[3][3],pom[3][3],g[3][3],z[3][3],rb[3][3],rp[3][3],rbp[3][3],p[3][3],q[3][3],qt[3][3];
    if(frame==TINYFRAMES_ITRF){tinyframes_erfa_Ir(out);return;}
    if(frame==TINYFRAMES_TEME||frame==TINYFRAMES_PEF1982) {
        tinyframes_erfa_Pom00(e->xp_rad,e->yp_rad,0,pom);
        if(frame==TINYFRAMES_PEF1982){copy(pom,out);return;}
        tinyframes_erfa_Ir(z);tinyframes_erfa_Rz(tinyframes_erfa_Gmst82(e->ut1[0],e->ut1[1]),z);multiply(pom,z,out);return;
    }
    model_matrices(e,model,ci,pom,g);
    if(frame==TINYFRAMES_GCRF){copy(g,out);return;}
    if(frame==TINYFRAMES_TIRS){copy(pom,out);return;}
    if(frame==TINYFRAMES_CIRS){tinyframes_erfa_Ir(z);tinyframes_erfa_Rz(tinyframes_erfa_Era00(e->ut1[0],e->ut1[1]),z);multiply(pom,z,out);return;}
    tinyframes_erfa_Bp00(e->tt[0],e->tt[1],rb,rp,rbp);
    if(frame==TINYFRAMES_J2000)copy(rb,q);
    else {
        if(frame==TINYFRAMES_MOD1976)tinyframes_erfa_Pmat76(e->tt[0],e->tt[1],p);
        else tinyframes_erfa_Pnm80(e->tt[0],e->tt[1],p);
        multiply(p,rb,q);
    }
    transpose(q,qt);multiply(g,qt,out);
}
static void matrix_unchecked(const tinyframes_epoch *e,tinyframes_frame from,tinyframes_frame to,tinyframes_model model,double out[3][3]) {
    double a[3][3],b[3][3],bt[3][3];
    if(from==to){tinyframes_erfa_Ir(out);return;}
    to_fixed(e,from,model,a);to_fixed(e,to,model,b);transpose(b,bt);multiply(bt,a,out);
}
tinyframes_status tinyframes_matrix(const tinyframes_epoch *epoch,tinyframes_frame from,tinyframes_frame to,tinyframes_model model,double out[3][3]) {
    tinyframes_epoch e;tinyframes_status status;double r[3][3];
    if(!out)return TINYFRAMES_INVALID;
    status=valid(epoch,from,to,model,&e);if(status!=TINYFRAMES_OK)return status;
    matrix_unchecked(&e,from,to,model,r);copy(r,out);return TINYFRAMES_OK;
}
tinyframes_status tinyframes_transform(const tinyframes_epoch *epoch,tinyframes_frame from,tinyframes_frame to,tinyframes_model model,const double r[3],const double v[3],double rout[3],double vout[3]) {
    tinyframes_epoch e,plus,minus;double m[3][3],mp[3][3],mm[3][3],rr[3]={0},vv[3]={0};tinyframes_status status;int i,j;
    if(!r||!v||!rout||!vout||rout==vout)return TINYFRAMES_INVALID;
    for(i=0;i<3;i++)if(!isfinite(r[i])||!isfinite(v[i]))return TINYFRAMES_INVALID;
    status=valid(epoch,from,to,model,&e);if(status!=TINYFRAMES_OK)return status;
    plus=e;minus=e;plus.tt[1]+=.5/86400.;minus.tt[1]-=.5/86400.;plus.ut1[1]+=.5*(1-e.lod_s/86400)/86400.;minus.ut1[1]-=.5*(1-e.lod_s/86400)/86400.;
    matrix_unchecked(&e,from,to,model,m);matrix_unchecked(&plus,from,to,model,mp);matrix_unchecked(&minus,from,to,model,mm);
    for(i=0;i<3;i++)for(j=0;j<3;j++){rr[i]+=m[i][j]*r[j];vv[i]+=m[i][j]*v[j]+(mp[i][j]-mm[i][j])*r[j];}
    for(i=0;i<3;i++)if(!isfinite(rr[i])||!isfinite(vv[i]))return TINYFRAMES_RANGE;
    for(i=0;i<3;i++){rout[i]=rr[i];vout[i]=vv[i];}return TINYFRAMES_OK;
}
