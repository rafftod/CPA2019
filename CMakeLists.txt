cmake_minimum_required(VERSION 3.0)

project(main)
find_package(Boost COMPONENTS program_options system thread REQUIRED)
include_directories(${BOOST_INCLUDE_DIR})

add_compile_options(-Wall -ggdb)

add_executable(main Index.cpp SequenceReader.cpp Header.cpp Smith_Waterman.cpp main.cpp)
TARGET_LINK_LIBRARIES(main LINK_PUBLIC ${Boost_LIBRARIES} Boost::thread )