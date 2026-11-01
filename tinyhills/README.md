# tinyhills

Heap-free C99 Clohessy–Wiltshire (Hill) linear relative-orbit propagation, 6×6 state transition matrices, and zero-terminal-position velocity targeting. Link `tinyhills::tinyhills`, include `tinyhills/tinyhills.h`.

Coordinates are x radial outward, y along-track, z chief-orbit normal. Positions are metres, relative rotating-frame velocities metres/second, mean motion radians/second, and time seconds. The chief orbit is circular and relative separation must be small compared with chief radius. This model is not a general nonlinear relative propagator or a relative-frame conversion library.

`tinyhills_propagate` advances a state; `tinyhills_transition` returns its row-major transition matrix ordered `[x,y,z,vx,vy,vz]`. `tinyhills_rendezvous` computes an initial velocity that reaches relative position zero. It rejects times where the terminal-position map is singular/nonunique, including the cross-track half-orbit degeneracies. `n=0` uses the constant-velocity limit. Valid profile: `n>=0`, finite times, `abs(n*dt)<=1e8`; all failures preserve outputs. Stack work is bounded and no mutable global state/heap is used.

Validation captures **30 propagation and 24 targeting cases** from Vallado Python `hillsr/hillsv`, pinned at `49df0479c950917c0a2c5ac31f2db6477bba57f6`. Three chief altitudes, varied signed initial states, negative/zero/positive intervals, time reversal and zero-terminal-position checks are covered. Maximum sampled differences on AppleClang 21 were **9.10e-13 m**, **2.23e-16 m/s**, and **6.89e-15 m/s** for targeting. Static tests also check singularity, invalid inputs, failure preservation, and n=0. Ordinary tests are offline; raw comparisons/generator are ignored under `.validation/tinyhills/`.

Implemented independently from the published CW transition equations; no Vallado source was translated. See LICENSE and NOTICE.
