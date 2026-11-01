cmake_minimum_required(VERSION 3.21)
project(tinyrigidbody VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyrigidbody SOURCES src/tinyrigidbody.c TEST_SOURCES tests/test_tinyrigidbody.c DEPENDS tinylinalg tinygyro)
