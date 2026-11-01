cmake_minimum_required(VERSION 3.21)
project(tinyhohmann VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyhohmann SOURCES src/tinyhohmann.c TEST_SOURCES tests/test_tinyhohmann.c DEPENDS tinyanomaly)
