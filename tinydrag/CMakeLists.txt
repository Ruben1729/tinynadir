cmake_minimum_required(VERSION 3.21)
project(tinydrag VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinydrag SOURCES src/tinydrag.c TEST_SOURCES tests/test_tinydrag.c)
