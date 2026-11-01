cmake_minimum_required(VERSION 3.21)
project(tinywheels VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinywheels SOURCES src/tinywheels.c TEST_SOURCES tests/test_tinywheels.c DEPENDS tinylinalg)
