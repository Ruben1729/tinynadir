cmake_minimum_required(VERSION 3.21)
project(tinycomplementary VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinycomplementary SOURCES src/tinycomplementary.c TEST_SOURCES tests/test_tinycomplementary.c DEPENDS tinylinalg tinygyro)
