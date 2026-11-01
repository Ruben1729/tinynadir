cmake_minimum_required(VERSION 3.21)
project(tinycovariance VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinycovariance SOURCES src/tinycovariance.c TEST_SOURCES tests/test_tinycovariance.c DEPENDS tinyelements)
