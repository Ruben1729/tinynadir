cmake_minimum_required(VERSION 3.21)
project(tinycoverage VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinycoverage SOURCES src/tinycoverage.c TEST_SOURCES tests/test_tinycoverage.c)
