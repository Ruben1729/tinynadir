cmake_minimum_required(VERSION 3.21)
project(tinyunscented VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyunscented SOURCES src/tinyunscented.c TEST_SOURCES tests/test_tinyunscented.c)
