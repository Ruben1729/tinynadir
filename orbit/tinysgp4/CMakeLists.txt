cmake_minimum_required(VERSION 3.21)
project(tinysgp4 VERSION 0.1.0 LANGUAGES C)

include(GNUInstallDirs)

set(_tgp4_native_default OFF)
if(PROJECT_IS_TOP_LEVEL AND NOT CMAKE_CROSSCOMPILING)
    set(_tgp4_native_default ON)
endif()
option(TGP4_BUILD_TESTS "Build offline reference-vector and API tests" ${_tgp4_native_default})
option(TGP4_BUILD_EXAMPLES "Build the prediction example" ${_tgp4_native_default})
option(TGP4_ENABLE_SANITIZERS "Enable AddressSanitizer and UBSan for native builds" OFF)
option(TGP4_INSTALL "Generate library/header/package installation rules" ${PROJECT_IS_TOP_LEVEL})

if(MSVC)
    set(_tgp4_math_default "")
else()
    set(_tgp4_math_default "m")
endif()
set(TGP4_MATH_LIBRARY "${_tgp4_math_default}" CACHE STRING
    "Math libraries to link; leave empty if supplied by the platform runtime")

if(TGP4_ENABLE_SANITIZERS)
    if(CMAKE_CROSSCOMPILING OR NOT CMAKE_C_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
        message(FATAL_ERROR "TGP4_ENABLE_SANITIZERS requires a native GCC or Clang build")
    endif()
endif()

function(tgp4_configure_target target)
    set_target_properties(${target} PROPERTIES
        C_STANDARD 99 C_STANDARD_REQUIRED YES C_EXTENSIONS NO)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4)
    elseif(CMAKE_C_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
        target_compile_options(${target} PRIVATE -Wall -Wextra)
    endif()
    if(TGP4_ENABLE_SANITIZERS)
        target_compile_options(${target} PRIVATE
            -fsanitize=address,undefined -fno-omit-frame-pointer)
    endif()
endfunction()

add_library(tinysgp4 STATIC src/sgp4.c src/tinysgp4.c)
add_library(tinysgp4::tinysgp4 ALIAS tinysgp4)
tgp4_configure_target(tinysgp4)
target_include_directories(tinysgp4 PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
if(TGP4_MATH_LIBRARY)
    target_link_libraries(tinysgp4 PUBLIC ${TGP4_MATH_LIBRARY})
endif()
if(TGP4_ENABLE_SANITIZERS)
    # Instrumented static objects require sanitizer runtimes in their consumers.
    target_link_options(tinysgp4 INTERFACE -fsanitize=address,undefined)
endif()

if(TGP4_BUILD_EXAMPLES)
    add_executable(tinysgp4_predict examples/predict.c)
    tgp4_configure_target(tinysgp4_predict)
    target_link_libraries(tinysgp4_predict PRIVATE tinysgp4::tinysgp4)
endif()

if(TGP4_BUILD_TESTS)
    enable_testing()
    add_executable(tinysgp4_verify tests/verify.c tests/data/cases.h
        tests/data/vallado.csv tests/data/published.csv)
    tgp4_configure_target(tinysgp4_verify)
    target_link_libraries(tinysgp4_verify PRIVATE tinysgp4::tinysgp4)
    add_test(NAME tinysgp4.reference COMMAND tinysgp4_verify)
    set_tests_properties(tinysgp4.reference PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
    if(TGP4_BUILD_EXAMPLES)
        add_test(NAME tinysgp4.example COMMAND tinysgp4_predict)
    endif()
    if(CMAKE_NM AND NOT MSVC)
        add_test(NAME tinysgp4.no_heap
            COMMAND "${CMAKE_COMMAND}"
                "-DNM=${CMAKE_NM}"
                "-DARCHIVE=$<TARGET_FILE:tinysgp4>"
                -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/CheckNoHeap.cmake")
    else()
        message(STATUS "Allocator-symbol test unavailable: no compatible nm tool")
    endif()
endif()

if(TGP4_INSTALL)
    include(CMakePackageConfigHelpers)
    set(_tgp4_package_dir "${CMAKE_INSTALL_LIBDIR}/cmake/tinysgp4")
    install(TARGETS tinysgp4 EXPORT tinysgp4Targets
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}")
    install(DIRECTORY include/tinysgp4 DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
    install(FILES LICENSE NOTICE README.md
        DESTINATION "${CMAKE_INSTALL_DATADIR}/doc/tinysgp4")
    install(EXPORT tinysgp4Targets FILE tinysgp4Targets.cmake
        NAMESPACE tinysgp4:: DESTINATION "${_tgp4_package_dir}")
    configure_package_config_file(cmake/tinysgp4Config.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/tinysgp4Config.cmake"
        INSTALL_DESTINATION "${_tgp4_package_dir}")
    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/tinysgp4ConfigVersion.cmake"
        VERSION "${PROJECT_VERSION}" COMPATIBILITY SameMajorVersion)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/tinysgp4Config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/tinysgp4ConfigVersion.cmake"
        DESTINATION "${_tgp4_package_dir}")
endif()
