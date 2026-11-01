cmake_minimum_required(VERSION 3.21)
project(tinysun VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinysun SOURCES
    src/tinysun.c
    src/erfa_epv00.c
    TEST_SOURCES tests/test_tinysun.c)
if(TINYSUN_INSTALL)
    install(FILES ERFA-LICENSE.txt DESTINATION "${CMAKE_INSTALL_DATADIR}/doc/tinysun")
endif()
