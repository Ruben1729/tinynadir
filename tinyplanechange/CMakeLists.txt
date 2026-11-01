cmake_minimum_required(VERSION 3.21)
project(tinyplanechange VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyplanechange SOURCES src/tinyplanechange.c TEST_SOURCES tests/test_tinyplanechange.c)
