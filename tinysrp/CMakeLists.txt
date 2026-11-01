cmake_minimum_required(VERSION 3.21)
project(tinysrp VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinysrp SOURCES src/tinysrp.c TEST_SOURCES tests/test_tinysrp.c)
