cmake_minimum_required(VERSION 3.21)
project(tinykepler VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinykepler SOURCES src/tinykepler.c TEST_SOURCES tests/test_tinykepler.c DEPENDS tinystumpff)
