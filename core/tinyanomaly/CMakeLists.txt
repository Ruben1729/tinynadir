cmake_minimum_required(VERSION 3.21)
project(tinyanomaly VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyanomaly SOURCES src/tinyanomaly.c TEST_SOURCES tests/test_tinyanomaly.c)
