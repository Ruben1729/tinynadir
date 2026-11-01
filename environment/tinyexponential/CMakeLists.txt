cmake_minimum_required(VERSION 3.21)
project(tinyexponential VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyexponential SOURCES src/tinyexponential.c TEST_SOURCES tests/test_tinyexponential.c)
