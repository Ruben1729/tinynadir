cmake_minimum_required(VERSION 3.21)
project(tinycowell VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinycowell SOURCES src/tinycowell.c TEST_SOURCES tests/test_tinycowell.c DEPENDS tinyrk4 tinyrkf45)
