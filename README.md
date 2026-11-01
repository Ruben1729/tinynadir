# TinyNadir

TinyNadir is a collection of algorithms for CubeSats and other embedded spacecraft systems. It brings established algorithms, such as Vallado's SGP4, into small C libraries with caller-owned storage and no heap allocation.

The collection is intended to grow to cover orbit propagation, attitude determination, and control. **The current implementation is the `tinysgp4` orbit propagation module.** Attitude estimators and controllers have not been implemented yet.

## Modules

| Module | Function | Status |
| --- | --- | --- |
| [tinysgp4](tinysgp4/README.md) | Double-precision SGP4, including near-Earth and deep-space propagation, from traditional TLEs | Implemented in C99; checked against frozen Vallado reference outputs |

Modules live in their own subdirectories with their source, headers, examples, tests, and attribution. They can be built independently or selected through the collection's CMake project. Test fixtures and examples are separate from the libraries used in firmware.

The design goals are explicit memory ownership, no hidden allocation, and documented numerical behavior and resource use. Memory and execution-time budgets depend on the chosen module, compiler, and platform. See each module's README for its supported inputs, validation coverage, and measured sizes.

## Build and test

Use CMake 3.21 or newer, a C99 compiler, and a supported build tool. From the repository root:

```sh
cmake -S . -B build/native -DCMAKE_BUILD_TYPE=Release
cmake --build build/native --config Release
ctest --test-dir build/native -C Release --output-on-failure
```

Native standalone builds enable module tests and examples by default. These default to disabled during cross-compilation or when TinyNadir is added to another project.

| Option | Default | Purpose |
| --- | --- | --- |
| `TINYNADIR_BUILD_TINYSGP4` | `ON` | Include the orbit propagation module. |
| `TINYNADIR_BUILD_TESTS` | Native standalone builds only | Enable tests for selected modules. |
| `TINYNADIR_BUILD_EXAMPLES` | Native standalone builds only | Enable examples for selected modules. |
| `TINYNADIR_INSTALL` | Standalone builds only | Generate installation rules for selected modules. |

Module options, such as `TGP4_BUILD_TESTS`, can override the collection defaults. For native AddressSanitizer and UBSan checks of the SGP4 module:

```sh
cmake -S . -B build/sanitize -DCMAKE_BUILD_TYPE=Debug -DTGP4_ENABLE_SANITIZERS=ON
cmake --build build/sanitize --config Debug
ctest --test-dir build/sanitize -C Debug --output-on-failure
```

To build just SGP4 independently, configure its directory:

```sh
cmake -S tinysgp4 -B build/sgp4 -DCMAKE_BUILD_TYPE=Release
cmake --build build/sgp4 --config Release
ctest --test-dir build/sgp4 -C Release --output-on-failure
```

## Use in firmware

Add the collection as a CMake subdirectory and link the modules your application needs:

```cmake
add_subdirectory(path/to/tinynadir)
target_link_libraries(app PRIVATE tinysgp4::tinysgp4)
```

The SGP4 public header remains `"tinysgp4/tinysgp4.h"`, and its API retains the `tgp4_` prefix. See the [module README](tinysgp4/README.md) and [prediction example](tinysgp4/examples/predict.c) for usage.

For cross-compilation, supply your platform toolchain file:

```sh
cmake -S . -B build/cross -DCMAKE_TOOLCHAIN_FILE=/path/to/toolchain.cmake \
    -DTINYNADIR_BUILD_TESTS=OFF -DTINYNADIR_BUILD_EXAMPLES=OFF
cmake --build build/cross
```

The toolchain selects the compiler, CPU/FPU flags, ABI, and platform runtime. The SGP4 module requires binary64 `double` and a math library; its README documents direct compilation without CMake.

## Installation

```sh
cmake --install build/native --config Release --prefix /path/to/prefix
```

Enabled modules install their own headers, libraries, CMake packages, and license notices. The current module can be consumed with `find_package(tinysgp4 CONFIG REQUIRED)` and the `tinysgp4::tinysgp4` target. Configure consumers with `-DCMAKE_PREFIX_PATH=/path/to/prefix`.

## License and attribution

TinyNadir contributions are available under the [MIT license](LICENSE). Individual modules preserve their upstream attribution and notices; see [tinysgp4/NOTICE](tinysgp4/NOTICE) for the SGP4 implementation and its original authors. Numerical comparisons, reference revisions, and the fixed-point investigation are documented in the [SGP4 README](tinysgp4/README.md).
