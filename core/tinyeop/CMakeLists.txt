cmake_minimum_required(VERSION 3.21)
project(tinyeop VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyeop SOURCES src/tinyeop.c TEST_SOURCES tests/test_tinyeop.c)
