cmake_minimum_required(VERSION 3.21)
project(tinymoon VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinymoon SOURCES
    src/tinymoon.c
    src/erfa_cp.c
    src/erfa_ir.c
    src/erfa_moon98.c
    src/erfa_obl06.c
    src/erfa_pfw06.c
    src/erfa_rx.c
    src/erfa_rxp.c
    src/erfa_rxpv.c
    src/erfa_rz.c
    src/erfa_s2pv.c
    TEST_SOURCES tests/test_tinymoon.c)
if(TINYMOON_INSTALL)
    install(FILES ERFA-LICENSE.txt DESTINATION "${CMAKE_INSTALL_DATADIR}/doc/tinymoon")
endif()
