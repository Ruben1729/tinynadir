cmake_minimum_required(VERSION 3.21)
project(tinyephemeris VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyephemeris SOURCES src/tinyephemeris.c TEST_SOURCES tests/test_tinyephemeris.c DEPENDS tinyinterpolation)
