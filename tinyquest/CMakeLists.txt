cmake_minimum_required(VERSION 3.21)
project(tinyquest VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyquest SOURCES src/tinyquest.c TEST_SOURCES tests/test_tinyquest.c DEPENDS tinylinalg tinyrotation)
