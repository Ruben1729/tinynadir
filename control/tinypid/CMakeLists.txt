cmake_minimum_required(VERSION 3.21)
project(tinypid VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinypid SOURCES src/tinypid.c TEST_SOURCES tests/test_tinypid.c DEPENDS tinypd)
