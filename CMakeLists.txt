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
option(TINYNADIR_ENABLE_SANITIZERS "Enable AddressSanitizer and UBSan in native module builds" OFF)

if(TINYNADIR_BUILD_TESTS)
    enable_testing()
endif()

if(TINYNADIR_BUILD_TINYSGP4)
    # Supply collection defaults without overriding explicit module options.
    # Normal variables let collection options take effect on reconfiguration too.
    foreach(_setting BUILD_TESTS BUILD_EXAMPLES INSTALL ENABLE_SANITIZERS)
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

# Modules are registered here after their independent comparison and offline
# regression tests pass. Additional modules are never mandatory dependencies
# of tinysgp4.
set(TINYNADIR_AVAILABLE_MODULES tinysgp4 tinylinalg tinystumpff tinyanomaly tinytime tinyrotation tinytriad tinydavenport tinyquest tinyigrf tinyeop tinyinterpolation tinyrk4 tinyrkf45 tinykepler tinyelements tinygeodesy tinythirdbody tinyexponential tinydrag tinysrp tinygravity)
foreach(_module IN LISTS TINYNADIR_AVAILABLE_MODULES)
    if(_module STREQUAL "tinysgp4")
        continue()
    endif()
    string(TOUPPER "${_module}" _upper)
    option(TINYNADIR_BUILD_${_upper} "Build the ${_module} module" ON)
    if(TINYNADIR_BUILD_${_upper} AND NOT TARGET ${_module}::${_module})
        add_subdirectory(${_module})
    endif()
endforeach()
