cmake_minimum_required(VERSION 3.21)
project(tinyhills VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyhills SOURCES src/tinyhills.c TEST_SOURCES tests/test_tinyhills.c)
