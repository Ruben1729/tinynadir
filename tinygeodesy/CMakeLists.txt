cmake_minimum_required(VERSION 3.21)
project(tinygeodesy VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinygeodesy SOURCES src/tinygeodesy.c TEST_SOURCES tests/test_tinygeodesy.c)
