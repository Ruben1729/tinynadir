cmake_minimum_required(VERSION 3.21)
project(tinygravity VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinygravity SOURCES src/tinygravity.c TEST_SOURCES tests/test_tinygravity.c)
