cmake_minimum_required(VERSION 3.21)
project(tinyinterpolation VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyinterpolation SOURCES src/tinyinterpolation.c TEST_SOURCES tests/test_tinyinterpolation.c)
