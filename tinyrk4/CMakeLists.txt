cmake_minimum_required(VERSION 3.21)
project(tinyrk4 VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyrk4 SOURCES src/tinyrk4.c TEST_SOURCES tests/test_tinyrk4.c)
