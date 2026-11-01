cmake_minimum_required(VERSION 3.21)
project(tinyrotation VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyrotation SOURCES src/tinyrotation.c TEST_SOURCES tests/test_tinyrotation.c)
