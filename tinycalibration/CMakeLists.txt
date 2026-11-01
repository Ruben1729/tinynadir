cmake_minimum_required(VERSION 3.21)
project(tinycalibration VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinycalibration SOURCES src/tinycalibration.c TEST_SOURCES tests/test_tinycalibration.c DEPENDS tinyrotation)
