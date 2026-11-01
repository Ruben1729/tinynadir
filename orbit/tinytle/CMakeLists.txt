cmake_minimum_required(VERSION 3.21)
project(tinytle VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinytle SOURCES src/tinytle.c src/initialize.c src/initialize_text.c TEST_SOURCES tests/test_tinytle.c DEPENDS tinytime tinysgp4)
