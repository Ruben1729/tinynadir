cmake_minimum_required(VERSION 3.21)
project(tinyspaceweather VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyspaceweather SOURCES src/tinyspaceweather.c TEST_SOURCES tests/test_tinyspaceweather.c)
