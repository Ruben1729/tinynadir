cmake_minimum_required(VERSION 3.21)
project(tinylaplace VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinylaplace SOURCES src/tinylaplace.c TEST_SOURCES tests/test_tinylaplace.c DEPENDS tinypolynomial tinylinalg)
