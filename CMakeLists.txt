cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 17)
project(automorphisms-bench)

add_subdirectory(lib/spqlios-arithmetic)

add_executable(main rlwe.cpp automorphisms.h automorphism.cpp main.cpp)
target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR}/lib/spqlios-arithmetic)
target_link_libraries(main libspqlios)
