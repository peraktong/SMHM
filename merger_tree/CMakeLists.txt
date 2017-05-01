cmake_minimum_required(VERSION 3.7)
project(merger_tree)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        erfinv.c
        gasdev.c
        get_redshift.c
        main.c
        main.h
        ran1.c
        zbrent.c)

add_executable(merger_tree ${SOURCE_FILES})