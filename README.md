# TinyNadir

TinyNadir provides heap-free C99 algorithms for CubeSat orbit prediction, environment modeling, attitude determination, and control. Libraries use double precision, caller-owned state/workspace, and bounded automatic storage.

**Start with a task, then select the pieces you need.** The six groups organize 68 algorithms; choosing one algorithm builds its required dependencies without pulling in its whole group.

## What are you trying to do?

| Task | Start here | A possible chain |
| --- | --- | --- |
| Detumble with magnetorquers | [Control](control/README.md) | Calibrated magnetic measurements → B-dot → torquer driver |
| Estimate attitude from Sun/magnetic measurements | [Attitude](attitude/README.md) | Common-frame reference vectors + calibrated body vectors → TRIAD or QUEST; optionally a temporal filter |
| Point with reaction wheels | [Control](control/README.md) | Attitude input + guidance → PD/PID → wheel allocation |
| Predict an orbit | [Orbit](orbit/README.md) | TLE/OMM → SGP4; choose Kepler or numerical propagation for suitable alternative initial conditions |
| Calculate a magnetic reference | [Environment](environment/README.md) | Position/time → IGRF → required reference frame |
| Fit an orbit, plan transfers or calculate access | [Mission](mission/README.md) | Observations → selected estimation/transfer/coverage helpers |

You do not need every piece. A B-dot detumbler needs no orbit or attitude estimate. A star tracker can replace vector-based attitude estimation. GNSS or uploaded states can replace an onboard propagator. The application supplies drivers, clocks, mission data, scheduling and mode handling.

## Browse by capability

| Group | Contains |
| --- | --- |
| [Core](core/README.md) | Math, rotations, time, coordinate systems, interpolation and integrators |
| [Orbit](orbit/README.md) | TLE/OMM inputs, SGP4, Kepler, J2 models and Cartesian propagation |
| [Environment](environment/README.md) | Sun, Moon, IGRF, gravity, atmospheres, forces and eclipses |
| [Attitude](attitude/README.md) | Calibration, TRIAD, QUEST, Davenport, gyro propagation and attitude filters |
| [Control](control/README.md) | Guidance, slews, B-dot, PD/PID, actuator allocation and rotational dynamics |
| [Mission](mission/README.md) | Orbit determination/fitting, uncertainty, transfers, rendezvous and coverage |

Each group README contains its algorithm glossary, selection advice and connections to other groups. Individual algorithm READMEs document APIs, supported profiles, memory requirements and numerical validation. These are explicit model profiles rather than a translation of every routine/variant in Vallado's repositories. SGP4-XP is outside this collection.

```text
tinynadir/
├── core/
├── orbit/          # e.g. orbit/tinysgp4/
├── environment/
├── attitude/       # e.g. attitude/tinyquest/
├── control/
├── mission/
├── examples/
└── cmake/
```

## Build and test

Use CMake 3.21+, a C99 compiler with eight-byte IEEE binary64 `double`, and your platform's math runtime. From the full checkout:

```sh
# All algorithms, offline tests and examples.
cmake -S . -B build/native -DCMAKE_BUILD_TYPE=Release
cmake --build build/native
ctest --test-dir build/native --output-on-failure

# Just a detumbler; no other group is included.
cmake -S . -B build/detumble -DCMAKE_BUILD_TYPE=Release \
  -DTINYNADIR_MODULES=tinybdot \
  -DTINYNADIR_BUILD_TESTS=OFF -DTINYNADIR_BUILD_EXAMPLES=OFF
cmake --build build/detumble
```

Select multiple algorithms with a quoted semicolon-separated list, such as `-DTINYNADIR_MODULES='tinypd;tinywheels'`. Required dependencies are found automatically across groups. Configure one algorithm directly with, for example, `cmake -S attitude/tinyquest -B build/quest`; keep the full checkout for shared helpers/dependencies. Existing build directories configured with the previous flat layout should be replaced with fresh build directories.

| Option | Default | Purpose |
| --- | --- | --- |
| `TINYNADIR_MODULES` | `all` | Algorithm names to build, plus required dependencies |
| `TINYNADIR_BUILD_TESTS` | Native top-level builds | Static regression, invariant and allocator-symbol checks |
| `TINYNADIR_BUILD_EXAMPLES` | Native top-level builds | Examples whose required algorithms are selected |
| `TINYNADIR_INSTALL` | Top-level builds | Installation rules for selected libraries/packages/notices |
| `TINYNADIR_ENABLE_SANITIZERS` | Off | Native GCC/Clang AddressSanitizer and UBSan |
| `TINYNADIR_BUILD_TINYNAME` | On | Omit an optional algorithm from an `all` build; dependencies may still require it |

Tests/examples default off for cross-compilation and subproject use. Supply your platform's `CMAKE_TOOLCHAIN_FILE` for firmware. Algorithm-specific test/install options remain available; SGP4 retains its historical `TGP4_*` options. For sanitizer checks, configure a separate Debug build with `-DTINYNADIR_ENABLE_SANITIZERS=ON`.

## Use in an application

Directory grouping does not change public headers, functions, targets or package names:

```cmake
set(TINYNADIR_MODULES "tinypd;tinywheels" CACHE STRING "Selected algorithms")
add_subdirectory(path/to/tinynadir EXCLUDE_FROM_ALL)
target_link_libraries(app PRIVATE tinypd::tinypd tinywheels::tinywheels)
```

Include headers such as `tinypd/tinypd.h`. SGP4 keeps `tinysgp4/tinysgp4.h` and its `tgp4_` API prefix. To use installed libraries:

```sh
cmake --install build/native --prefix /path/to/prefix
```

Set the consumer's `CMAKE_PREFIX_PATH`, call `find_package(tinyquest CONFIG REQUIRED)`, and link `tinyquest::tinyquest`. Required packages are found transitively. There is no aggregate group library or monolithic `find_package(tinynadir)` requirement.

Runnable examples demonstrate [minimal detumbling](examples/detumble.c), [wheel pointing from an attitude input](examples/wheel_pointing.c), and [SGP4/Sun/IGRF/TRIAD snapshot attitude](examples/snapshot_attitude.c). The snapshot uses illustrative orbit/EOP data and synthetic sensors to demonstrate connections and conventions.

## Conventions, footprint and validation

Most new state/force interfaces use SI and radians. Legacy SGP4 outputs are **km/km/s**, propagation time is **minutes**, and GP input/fitting fields retain **degrees, revolutions/day and B* in inverse Earth radii**. Attitude quaternions are **scalar-first Hamilton, reference-to-body**. Check units, frames, epochs and model definitions at every connection; TLE mean elements are not ordinary osculating Kepler elements.

There is no hidden clock, EOP/leap database, space-weather service or environment archive. Supply applicable data and honor model ranges. Workspace figures exclude additional stack, callbacks and platform runtime; measure the linked selection on your target. Group READMEs include representative resource figures, and [SGP4 documents the evaluated fixed-point approach](orbit/tinysgp4/README.md#fixed-point-investigation) that was dropped in favor of binary64.

All 68 algorithms have static offline regression fixtures and relevant invariant/failure checks. Comparisons used pinned Vallado implementations and independent equivalents such as SciPy, ERFA, AHRS, FilterPy, ppigrf, original NRLMSISE-00 and Orekit DTM2000. Exact versions/commits, tolerances and sampled differences remain in the algorithm READMEs/notices, including the [SGP4 provenance table](orbit/tinysgp4/README.md#compared-implementations-and-data-provenance). Reference code/packages/generators/raw results remain ignored under `.validation/`; tests need no external algorithm library or Python/Java environment.

A Cortex-M7 compile/archive check with Arm GCC 16.2, `-Os`, FPv5-D16 hardware double precision and genuine Newlib headers built all 68 algorithms as strict C99. All archives had zero `.data`/`.bss` and no direct allocator imports. Unlinked code plus read-only data totaled **342,160 bytes** for every algorithm together, excluding libm/libc/libgcc, startup, application code, caller storage, stack and linker garbage collection. This is compile/storage evidence; hardware timing and final firmware size require target measurements.

Reference agreement verifies selected implementations/cases. It does not establish physical prediction accuracy, sensor performance, closed-loop stability or flight qualification.

## License

Original TinyNadir contributions use [MIT](LICENSE). Selected algorithms retain additional upstream terms: ERFA BSD-style notices, DTM2000 Apache-2.0, NRLMSISE-00 public-domain provenance, and IGRF coefficient CC BY 4.0 attribution. Read the LICENSE/NOTICE files for the algorithms you use; installation preserves those notices.
