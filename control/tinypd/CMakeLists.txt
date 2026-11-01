cmake_minimum_required(VERSION 3.21)
project(tinypd VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinypd SOURCES src/tinypd.c TEST_SOURCES tests/test_tinypd.c DEPENDS tinyrotation)
