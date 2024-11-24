cmake_minimum_required(VERSION 3.20)
project(Iterative-methods)

enable_testing()

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
    set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -pedantic")
endif()

set(SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/src/matrix.cpp)
add_executable(main main.cpp ${SOURCE})

set(INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(main PRIVATE ${INCLUDE})

option(ERRORS "display errors" OFF)
if(ERRORS)
    target_compile_definitions(main PRIVATE DISPLAY_ERRORS)
endif()