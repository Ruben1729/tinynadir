cmake_minimum_required(VERSION 3.21)
project(tinytriad VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinytriad SOURCES src/tinytriad.c TEST_SOURCES tests/test_tinytriad.c DEPENDS tinylinalg tinyrotation)
