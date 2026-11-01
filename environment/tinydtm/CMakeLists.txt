cmake_minimum_required(VERSION 3.21)
project(tinydtm VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinydtm SOURCES src/tinydtm.c TEST_SOURCES tests/test_tinydtm.c)
