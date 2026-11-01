cmake_minimum_required(VERSION 3.21)
project(tinygauss VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinygauss SOURCES src/tinygauss.c TEST_SOURCES tests/test_tinygauss.c DEPENDS tinypolynomial tinylinalg tinygibbs)
