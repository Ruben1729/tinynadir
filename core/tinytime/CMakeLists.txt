cmake_minimum_required(VERSION 3.21)
project(tinytime VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinytime SOURCES src/tinytime.c TEST_SOURCES tests/test_tinytime.c)
