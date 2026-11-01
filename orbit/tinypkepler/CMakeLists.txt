cmake_minimum_required(VERSION 3.21)
project(tinypkepler VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinypkepler SOURCES src/tinypkepler.c TEST_SOURCES tests/test_tinypkepler.c DEPENDS tinyelements tinyanomaly)
