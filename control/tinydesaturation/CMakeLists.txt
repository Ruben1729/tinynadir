cmake_minimum_required(VERSION 3.21)
project(tinydesaturation VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinydesaturation SOURCES src/tinydesaturation.c TEST_SOURCES tests/test_tinydesaturation.c DEPENDS tinymagnetic)
