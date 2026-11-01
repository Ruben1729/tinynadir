# tinyshadow

Spherical-body sunlight fraction from overlap of apparent angular discs, cylindrical parallel-light shadow, ellipsoid segment visibility, and a bounded bracketed threshold-event solver. Sun fraction includes full light, umbra, penumbra and annular geometry in the stated angular-disc approximation; it is not an atmospheric/refraction or exact surface-brightness model. Positions/axes use SI in one common frame/origin. Event detection finds one already-bracketed crossing; the application supplies scanning and avoids missing short events. No heap/global state; only bounded local vectors/scalars. Outputs preserved on failure.

## Verification

120 independent SciPy adaptive numerical integrals of apparent-disc intersections validate fraction to5e-8 absolute;128 SciPy bounded-minimization ellipsoid segment cases validate visibility. A SciPy Brent threshold event and analytic cylindrical cases add independent checks. Raw numerical integration tolerances/cases and maximum errors stay in ignored validation. References: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html and https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html .

Reference scripts, packages and raw outputs are in ignored `.validation/tinyshadow/`. Regression fixtures are checked in and run offline without Python or reference libraries.
