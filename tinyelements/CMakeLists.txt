cmake_minimum_required(VERSION 3.21)
project(tinyelements VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyelements SOURCES src/tinyelements.c TEST_SOURCES tests/test_tinyelements.c)
