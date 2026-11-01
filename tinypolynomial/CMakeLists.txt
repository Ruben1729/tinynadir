cmake_minimum_required(VERSION 3.21)
project(tinypolynomial VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinypolynomial SOURCES src/tinypolynomial.c TEST_SOURCES tests/test_tinypolynomial.c)
