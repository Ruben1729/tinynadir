cmake_minimum_required(VERSION 3.21)
project(tinydoubler VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinydoubler SOURCES src/tinydoubler.c TEST_SOURCES tests/test_tinydoubler.c DEPENDS tinylinalg tinygibbs)
