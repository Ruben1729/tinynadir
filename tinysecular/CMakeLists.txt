cmake_minimum_required(VERSION 3.21)
project(tinysecular VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinysecular SOURCES src/tinysecular.c TEST_SOURCES tests/test_tinysecular.c)
