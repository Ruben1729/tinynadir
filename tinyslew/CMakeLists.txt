cmake_minimum_required(VERSION 3.21)
project(tinyslew VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyslew SOURCES src/tinyslew.c TEST_SOURCES tests/test_tinyslew.c DEPENDS tinyrotation)
