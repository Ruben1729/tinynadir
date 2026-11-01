# tinysgp4

Part of [Orbit](../README.md) in [TinyNadir](../../README.md).

tinysgp4 is a heap-free C99 library for predicting satellite position and velocity from a two-line element set (TLE). It implements near-Earth and deep-space SGP4 using double precision.

It is the orbit propagation module of [TinyNadir](../../README.md). This directory contains the standalone library, examples, tests, and upstream notices. The public C API and CMake target retain their original names.

The library uses caller-owned storage, has no mutable global state, and performs no file, console, or network I/O. Building and testing require no downloaded source, external SGP4 library, C++ compiler, or Python. Library source and test data are checked in; CMake generates only build files. The standard C math library is required.

## Build

Use **CMake 3.21 or newer**, a C99 compiler, and a supported build tool such as Ninja. CMake builds the library, examples, and tests in C99 mode with language extensions disabled. All configuration and tests run offline.

Run the commands in this README from the `tinysgp4/` module directory (`cd orbit/tinysgp4` from the collection root). To build through the collection, follow the [root build instructions](../../README.md#build-and-test).

```sh
cmake -S . -B build/native -DCMAKE_BUILD_TYPE=Release
cmake --build build/native --config Release
ctest --test-dir build/native -C Release --output-on-failure
```

This builds the static library, the `tinysgp4_predict` example, and the `tinysgp4_verify` test executable. CTest also checks the library for direct allocator dependencies when a compatible `nm` tool is available. Tests read their checked-in vectors using an explicit source-directory working path, so they can be run from any build directory.

Tests and examples default to enabled for standalone native builds, and disabled for cross-compilation or use through `add_subdirectory`. Installation rules default to enabled only for standalone builds. Configure these options with `-D`:

| Option | Purpose |
| --- | --- |
| `TGP4_BUILD_TESTS` | Build the offline reference-vector and API tests. |
| `TGP4_BUILD_EXAMPLES` | Build the prediction example. |
| `TGP4_ENABLE_SANITIZERS` | Enable AddressSanitizer and UBSan on native GCC/Clang builds; off by default. |
| `TGP4_INSTALL` | Generate header, library, and CMake package installation rules. |
| `TGP4_MATH_LIBRARY` | Platform math libraries; normally `m`, or empty for MSVC. Set empty if your runtime already supplies them. |

### Use from another project

For a source checkout, link the exported build target:

```cmake
add_subdirectory(path/to/tinynadir/orbit/tinysgp4)
target_link_libraries(app PRIVATE tinysgp4::tinysgp4)
```

Alternatively, install the library and headers:

```sh
cmake --install build/native --config Release --prefix /path/to/prefix
```

Configure a consumer with `-DCMAKE_PREFIX_PATH=/path/to/prefix`, then use:

```cmake
find_package(tinysgp4 CONFIG REQUIRED)
target_link_libraries(app PRIVATE tinysgp4::tinysgp4)
```

The target supplies its include path and math-library link dependency. Consumer projects do not need to set these manually. Installation includes `LICENSE` and `NOTICE`.

### Sanitizers and optimization checks

```sh
cmake -S . -B build/sanitize -DCMAKE_BUILD_TYPE=Debug -DTGP4_ENABLE_SANITIZERS=ON
cmake --build build/sanitize --config Debug
ctest --test-dir build/sanitize -C Debug --output-on-failure
```

For a GCC/Clang optimization matrix, use separate build directories:

```sh
for level in 0 2 3; do
    cmake -S . -B "build/O${level}" -DCMAKE_BUILD_TYPE=Release \
        "-DCMAKE_C_FLAGS_RELEASE=-O${level} -DNDEBUG"
    cmake --build "build/O${level}" --config Release
    ctest --test-dir "build/O${level}" -C Release --output-on-failure
done
```

### Cross-compilation or direct compilation

Supply compiler and platform settings in your own CMake toolchain file:

```sh
cmake -S . -B build/cross -DCMAKE_TOOLCHAIN_FILE=/path/to/toolchain.cmake \
    -DTGP4_BUILD_TESTS=OFF -DTGP4_BUILD_EXAMPLES=OFF
cmake --build build/cross
```

The toolchain chooses the compiler, CPU/FPU flags, ABI, and platform libraries. Bare-metal toolchains can set `CMAKE_TRY_COMPILE_TARGET_TYPE` to `STATIC_LIBRARY` when compiler checks cannot link executables without board startup code. Configure each toolchain in its own build directory. Enable floating-point hardware in platform startup where applicable, and use compatible floating-point calling conventions throughout your application.

To integrate without CMake, compile `src/sgp4.c` and `src/tinysgp4.c` in C99 mode (`-std=c99` on GCC/Clang), add `include/` to the header search path, and link your platform's math library (`-lm` on typical Unix-like systems). The implementation requires an eight-byte binary64 `double`; this is checked at compile time.

## Usage

```c
#include "tinysgp4/tinysgp4.h"

const char line1[] =
    "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
const char line2[] =
    "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";

tgp4_satellite satellite = {0};
tgp4_state state;
tgp4_status status = tgp4_initialize(
    &satellite, line1, sizeof line1 - 1, line2, sizeof line2 - 1);

if (status == TGP4_OK) {
    status = tgp4_propagate(&satellite, 60.0, &state);
    if (status == TGP4_OK) {
        /* state.position_km[3] and state.velocity_km_s[3] */
    }
}
```

Initialize once per TLE, then reuse the satellite for subsequent requests. The time argument is **minutes since the TLE epoch**; negative values request earlier states. Outputs are Cartesian vectors in **TEME** (True Equator, Mean Equinox), with position in **kilometres** and velocity in **kilometres per second**. Ground coordinates, antenna look angles, and absolute-date conversion belong to the application.

The API currently accepts finite times within **±10080 minutes (seven days)**. This is an interface limit, not a guarantee that a TLE remains accurate for seven days. Prediction quality depends on TLE age, maneuvers, drag, and the application.

### Input and state rules

- Supply the two 69-character element lines, excluding the optional satellite name. Buffers need not be NUL-terminated; trailing CR/LF is accepted. Pass the actual length, excluding a terminating NUL.
- Both lines must have valid checksums and matching five-digit numeric satellite IDs. Alpha-5 identifiers and SGP4-XP element sets are not supported.
- The library uses **WGS-72** gravity constants and improved SGP4 operation mode (`i`).
- TLE years `00`–`56` map to 2000–2056; `57`–`99` map to 1957–1999. Leap-second conversion is outside the API.
- Zero-initialize `tgp4_satellite` before first use. Failed initialization or propagation preserves the satellite; failed propagation also preserves the output.
- Check the returned status before reading an output. Different satellite instances are independent; serialize calls that access the same instance.
- Fields exposed through `internal.h` are implementation details. Use the public functions to manage a satellite.

See [examples/predict.c](examples/predict.c) for a complete program.

## Resource use

These measurements describe the library, excluding the test corpus and external math/runtime implementations.

| Resource | Measured size |
| --- | ---: |
| One satellite instance, 32-bit ABI | 1,560 bytes |
| One satellite instance, tested 64-bit ABI | 1,568 bytes |
| One output state | 48 bytes |
| Required heap allocation | 0 bytes |
| Mutable library globals | 0 bytes |
| Library stack during initialization, measured call chain | 2,768 bytes |
| Library stack during propagation, measured call chain | 2,272 bytes |
| Library code and read-only constants | 17,024 bytes, about 16.6 KiB |

The code-size and stack measurements used a 32-bit Arm ABI, GCC 16.2, `-O2`, and hardware double precision. Sizes vary with ABI, compiler, optimization, and platform libraries; use `sizeof(tgp4_satellite)` for your build's per-satellite storage.

Stack figures sum compiler-reported library frames. They exclude external math/runtime frames, caller variables, interrupts, and operating-system context. Add satellite/output storage if those objects are local variables. Budget approximately **1.5 KiB per stored satellite**, plus working stack and your platform's overhead. The library objects have no initialized or zero-initialized mutable global storage.

The final linked image can add math/runtime code or remove unused functions. Measure final flash use, stack high-water use, and execution time in your application. No universal timing or total-stack bound is claimed.

## Verification

The reference comparisons were performed during development. Their expected inputs, states, and errors are now checked in under `tests/data/`. Normal tests replay these static fixtures entirely offline; they do not compile or fetch another implementation.

| Checked-in coverage | Count |
| --- | ---: |
| TLE fixtures, including an intentional initialization failure and repeated satellite | 34 |
| Successful states captured from Vallado's implementation | 19,530 |
| Captured requests expected to report a propagation error | 1,400 |
| Additional published verification states in legacy mode | 597 |

The 1,400 expected errors are **passing tests when the library rejects those requests**. They are not outstanding test failures. One fixture intentionally fails initialization. Test vectors are used only by the test executable and are not linked into the library or application.

The host tests read the local CSV files and require no reference executable. The recorded sequence covers a 481-point grid over ±7 days, probes around ±720 minutes, negative and reversed requests, repeated times, deterministic randomized requests, and interleaved satellites. It therefore continues to exercise mutable-cache behavior without a running reference library. Converted TLE input values are recorded as additional regression checks.

Successful states use absolute tolerances of **1 mm in position** and **0.00001 m/s in velocity**, plus **32 binary64 epsilons times the expected vector norm**. That small relative allowance accommodates rounding at the enormous coordinates in some synthetic stress cases; for ordinary orbital distances the absolute tolerance dominates. In the development comparison with Vallado, the largest observed differences were approximately **0.14 micrometres** and **8.43 × 10⁻¹¹ m/s**. These are sampled implementation differences, not satellite-tracking accuracy guarantees. The older published vectors are less precise; their test explicitly selects legacy mode (`a`) while the public API remains in improved mode.

API regressions cover truncated and nonterminated buffers, checksums, invalid fields, null pointers, nonfinite times, failure-state preservation, the TLE century pivot, leap years, midnight, noon, and day boundaries. The exact TLE from [vaelen/sgp4 PR #1](https://github.com/vaelen/sgp4/pull/1), opened September 1, 2026, verifies preservation of its epoch. That report concerned a half-day shift in TLE formatting; this library has no TLE formatter.

### Compared implementations and data provenance

| Repository | Pinned commit | Purpose |
| --- | --- | --- |
| [CelesTrak/fundamentals-of-astrodynamics](https://github.com/CelesTrak/fundamentals-of-astrodynamics/tree/49df0479c950917c0a2c5ac31f2db6477bba57f6) | `49df0479c950917c0a2c5ac31f2db6477bba57f6` | Direct comparison with David Vallado's own C++ SGP4, version **2025-12-30**, using its TLE parser and propagator. The source is linked by [Vallado's software page](https://celestrak.org/software/vallado-sw.php). |
| [aholinch/sgp4](https://github.com/aholinch/sgp4/tree/552cb1489a52c3023ae70cb6c7e239e84c5950fe) | `552cb1489a52c3023ae70cb6c7e239e84c5950fe` | Original C adaptation used as the implementation starting point and initial comparison reference; baseline **2016-03-09**. |

`tests/data/vallado.csv` contains recorded outputs of the official C++ `twoline2rv` and `sgp4`, using WGS-72 and improved mode. It preserves all 20,930 sampled requests, including the 1,400 error outcomes. `tests/data/cases.h` records their normalized TLE inputs and converted elements. Only invalid verification-file checksums and blank ephemeris-type fields were normalized before numerical comparison.

`tests/data/published.csv` retains the 597 applicable rows from Vallado's `tcppver.out`, distributed in the pinned aholinch repository. The source corpus was `SGP4-VER.TLE` from the [historical CelesTrak distribution](https://celestrak.org/publications/AIAA/2006-6753/). Of the original published output rows, 69 exceeded this API's time range and one belonged to the intentional initialization failure.

The direct comparison exposed amplified rounding differences in TLE unit conversion. The parser now follows Vallado's evaluation order for degrees, mean motion, and B*. This removed a 73.94 m synthetic-case outlier without changing the orbital equations or relaxing test tolerances. Expected values were captured from the reference implementations, not calculated from tinysgp4 itself. No reference code, download machinery, or C++ bridge is required for the retained tests.

## Fixed-point investigation

A signed 64-bit Q19.44 fixed-point implementation was evaluated and set aside:

| Sampled interval | Worst position difference from the double reference |
| --- | ---: |
| At the TLE epoch | 1.14 m |
| Within ±1 day | 426 m |
| Entire sample through +7 days | 12.44 km |

The fixed scale lost small coefficients and also produced an extra arithmetic failure. These results did not justify keeping that backend. **The library uses double precision; the abandoned fixed-point implementation has been removed.** This conclusion concerns the evaluated implementation, not every possible fixed-point design.

## License and attribution

[MIT](LICENSE). Commercial use, modification, redistribution, and sublicensing are permitted; retain the license and copyright notice. The original upstream attribution and Unlicense notice are preserved in [NOTICE](NOTICE).

Credit: David A. Vallado, Paul Crawford, Richard Hujsak, and T. S. Kelso, [*Revisiting Spacetrack Report #3*, AIAA 2006-6753](https://celestrak.org/publications/AIAA/2006-6753/), and aholinch for the original C adaptation.
