cmake_minimum_required(VERSION 3.21)
project(tinylambert VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinylambert SOURCES src/tinylambert.c TEST_SOURCES tests/test_tinylambert.c DEPENDS tinystumpff)
