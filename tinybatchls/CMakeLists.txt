cmake_minimum_required(VERSION 3.21)
project(tinybatchls VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinybatchls SOURCES src/tinybatchls.c TEST_SOURCES tests/test_tinybatchls.c)
