cmake_minimum_required(VERSION 3.21)
project(tinygibbs VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinygibbs SOURCES src/tinygibbs.c TEST_SOURCES tests/test_tinygibbs.c)
