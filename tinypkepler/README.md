# tinypkepler

Heap-free C99 **Vallado-style first-order J2 secular-precession profile with ndot=nddot=0**. Include `tinypkepler/tinypkepler.h`; link `tinypkepler::tinypkepler` (uses `tinyelements` and `tinyanomaly`).

Inputs are mu (m³/s²), equatorial radius (m), dimensionless J2, model-specific mean classical elements `[p,e,i,Omega,omega,nu]`, and elapsed seconds. The profile keeps p/e/i constant, advances node at `-1.5*n*J2*(Re/p)^2*cos(i)`, periapsis at `.75*n*J2*(Re/p)^2*(5*cos(i)^2-1)`, and mean anomaly at Keplerian `n=sqrt(mu/a^3)`. It deliberately matches this pkepler profile's mean-rate convention; it does **not** add a separate J2 mean-anomaly-rate correction.

Only elliptic e<1 is supported; `abs(n*dt)<=1e8`. The initial elements belong to this approximate secular theory, not SGP4/TLE mean elements. There is no J4, drag/mean-motion derivative, or short-period reconstruction in this implementation. Selecting this approximation does not improve a TLE or replace a force-model numerical propagator. J2=0 gives the corresponding two-body mean-element evolution.

Output includes the evolved model elements, reconstructed inertial state, and secular rates. No allocation or mutable globals are used; failure preserves output.

Validation captures **36 pinned Vallado Python pkepler cases**, commit `49df0479c950917c0a2c5ac31f2db6477bba57f6`, with ndot=nddot=0, eccentricities .01/.2/.6, inclinations .2/1.1/2.4, and times -1 day/zero/1 hour/+1 day. Maximum sampled component differences on AppleClang 21 were **3.06e-7 m** and **1.86e-10 m/s**. Static limits are 1 mm and 1e-6 m/s; constant p/e/i, J2=0, rejected conics, and preservation are checked. Oracle source/raw output/generator stay ignored under `.validation/tinypkepler/`; normal tests are offline.

Independent implementation from published secular precession and Kepler equations, consistent with the documented Algorithm66 profile. No AGPL source was translated. See LICENSE/NOTICE.
