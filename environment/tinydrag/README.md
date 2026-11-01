# tinydrag

Part of [Environment](../README.md) in [TinyNadir](../../README.md).

Cannonball aerodynamic drag acceleration from explicit atmospheric density, spacecraft-relative air velocity, Cd, area and mass in SI units. The caller supplies atmospheric rotation/winds and velocity frame conversion; no density model or SGP4 dependency is imposed. No lift, surface-panel aerodynamics or implicit B* conversion. No heap, mutable globals or I/O; three-double output scratch plus scalars.

## Verification

100 independent NumPy vector-norm force evaluations are frozen, with componentwise relative tolerance 5e-15 and a drag-work sign check. Covers varied densities, areas/masses, directions and zero density. Reference formula and existing equivalent force model description: https://help.agi.com/stkEngineOnUNIX/Content/hpop/hpop.htm ; numerical oracle https://numpy.org/doc/stable/reference/generated/numpy.linalg.norm.html .

Reference scripts, packages and raw outputs are in ignored `.validation/tinydrag/`. Regression fixtures are checked in and run offline without Python or reference libraries.
