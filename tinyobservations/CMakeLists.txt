cmake_minimum_required(VERSION 3.21)
project(tinyobservations VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyobservations SOURCES src/tinyobservations.c TEST_SOURCES tests/test_tinyobservations.c DEPENDS tinygeodesy)
