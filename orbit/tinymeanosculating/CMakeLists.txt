cmake_minimum_required(VERSION 3.21)
project(tinymeanosculating VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinymeanosculating SOURCES src/tinymeanosculating.c TEST_SOURCES tests/test_tinymeanosculating.c DEPENDS tinycowell tinyelements tinyanomaly)
