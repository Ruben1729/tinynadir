cmake_minimum_required(VERSION 3.21)
project(tinyomm VERSION 0.1.0 LANGUAGES C)
include(../../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyomm SOURCES src/tinyomm.c src/initialize.c TEST_SOURCES tests/test_tinyomm.c DEPENDS tinytle tinytime)
