cmake_minimum_required(VERSION 3.21)
project(tinyrendezvous VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinyrendezvous SOURCES src/tinyrendezvous.c TEST_SOURCES tests/test_tinyrendezvous.c DEPENDS tinykepler tinylambert)
