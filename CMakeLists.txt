cmake_minimum_required(VERSION 3.0)


set(Boost_USE_STATIC_LIBS   ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_ADDITIONAL_VERSIONS "1.44" "1.44.0")


project(main)
find_package(Threads REQUIRED)



add_compile_options(-Wall -ggdb)

add_executable(main Index.cpp SequenceReader.cpp Header.cpp Smith_Waterman.cpp main.cpp)
