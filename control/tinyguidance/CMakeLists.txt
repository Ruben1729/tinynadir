cmake_minimum_required(VERSION 3.21)
project(tinyguidance VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyguidance SOURCES src/tinyguidance.c TEST_SOURCES tests/test_tinyguidance.c DEPENDS tinytriad tinylinalg)
