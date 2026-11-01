cmake_minimum_required(VERSION 3.21)
project(tinynadir VERSION 0.1.0 LANGUAGES C)

set(_tinynadir_native_default OFF)
if(PROJECT_IS_TOP_LEVEL AND NOT CMAKE_CROSSCOMPILING)
    set(_tinynadir_native_default ON)
endif()

option(TINYNADIR_BUILD_TINYSGP4 "Build the tinysgp4 orbit propagation module" ON)
option(TINYNADIR_BUILD_TESTS "Build tests for enabled modules" ${_tinynadir_native_default})
option(TINYNADIR_BUILD_EXAMPLES "Build examples for enabled modules" ${_tinynadir_native_default})
option(TINYNADIR_INSTALL "Generate installation rules for enabled modules" ${PROJECT_IS_TOP_LEVEL})

if(TINYNADIR_BUILD_TINYSGP4)
    # Supply collection defaults without overriding explicit module options.
    # Normal variables let collection options take effect on reconfiguration too.
    foreach(_setting BUILD_TESTS BUILD_EXAMPLES INSTALL)
        if(NOT DEFINED TGP4_${_setting})
            set(TGP4_${_setting} "${TINYNADIR_${_setting}}")
        endif()
    endforeach()

    # CTest needs to be enabled here to discover tests in the module directory.
    if(TGP4_BUILD_TESTS)
        enable_testing()
    endif()
    add_subdirectory(tinysgp4)
endif()
