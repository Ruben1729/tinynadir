cmake_minimum_required(VERSION 3.21)
project(tinysgp4fit VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinysgp4fit SOURCES src/tinysgp4fit.c
    TEST_SOURCES tests/test_tinysgp4fit.c DEPENDS tinybatchls tinytle)
