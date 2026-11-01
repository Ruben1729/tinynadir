cmake_minimum_required(VERSION 3.21)
project(tinymekf VERSION 0.1.0 LANGUAGES C)
include(../cmake/TinyNadirModule.cmake)
tinynadir_add_module(NAME tinymekf SOURCES src/tinymekf.c TEST_SOURCES tests/test_tinymekf.c DEPENDS tinylinalg tinygyro)
