cmake_minimum_required(VERSION 3.21)
project(tinyigrf VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyigrf SOURCES src/tinyigrf.c TEST_SOURCES tests/test_tinyigrf.c)
