cmake_minimum_required(VERSION 3.0)

project(main)

add_compile_options(-Wall -ggdb)

add_executable(main Index.cpp SequenceReader.cpp Header.cpp main.cpp)
