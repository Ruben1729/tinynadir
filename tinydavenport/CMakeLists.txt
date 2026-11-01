cmake_minimum_required(VERSION 3.21)
project(tinydavenport VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinydavenport SOURCES src/tinydavenport.c TEST_SOURCES tests/test_tinydavenport.c DEPENDS tinylinalg tinyrotation)
