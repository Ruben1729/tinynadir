cmake_minimum_required(VERSION 3.21)
project(tinyrkf45 VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyrkf45 SOURCES src/tinyrkf45.c TEST_SOURCES tests/test_tinyrkf45.c)
