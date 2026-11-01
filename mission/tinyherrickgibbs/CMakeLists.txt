cmake_minimum_required(VERSION 3.21)
project(tinyherrickgibbs VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyherrickgibbs SOURCES src/tinyherrickgibbs.c TEST_SOURCES tests/test_tinyherrickgibbs.c)
