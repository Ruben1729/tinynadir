cmake_minimum_required(VERSION 3.21)
project(tinymagnetic VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinymagnetic SOURCES src/tinymagnetic.c TEST_SOURCES tests/test_tinymagnetic.c DEPENDS tinylinalg)
