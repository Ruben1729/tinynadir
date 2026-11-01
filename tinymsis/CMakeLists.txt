cmake_minimum_required(VERSION 3.21)
project(tinymsis VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinymsis SOURCES src/tinymsis.c src/core.c src/coefficients.c TEST_SOURCES tests/test_tinymsis.c)
