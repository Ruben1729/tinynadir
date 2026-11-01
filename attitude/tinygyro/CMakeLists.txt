cmake_minimum_required(VERSION 3.21)
project(tinygyro VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinygyro SOURCES src/tinygyro.c TEST_SOURCES tests/test_tinygyro.c DEPENDS tinyrotation)
