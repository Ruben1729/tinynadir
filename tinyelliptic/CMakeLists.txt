cmake_minimum_required(VERSION 3.21)
project(tinyelliptic VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyelliptic SOURCES src/tinyelliptic.c TEST_SOURCES tests/test_tinyelliptic.c)
