cmake_minimum_required(VERSION 3.20)
project(epsilon_active_placer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(OpenMP REQUIRED)

add_library(epsilon_active_core
    src/batch_acceptance.cpp
    src/bisection.cpp
    src/bookshelf.cpp
    src/coarse_flow.cpp
    src/compact_recovery.cpp
    src/density.cpp
    src/hpwl.cpp
    src/lambda_controller.cpp
    src/optimizer.cpp
    src/placer.cpp
    src/recovery.cpp
    src/swap_recovery.cpp
    src/transport.cpp)
target_include_directories(epsilon_active_core PUBLIC include)
target_link_libraries(epsilon_active_core PUBLIC OpenMP::OpenMP_CXX)

if(MSVC)
    target_compile_options(epsilon_active_core PRIVATE /W4 /O2)
else()
    target_compile_options(epsilon_active_core PRIVATE -Wall -Wextra -Wpedantic -O3 -march=native)
endif()

add_executable(epsilon_active src/main.cpp)
target_link_libraries(epsilon_active PRIVATE epsilon_active_core)

enable_testing()
add_executable(core_tests tests/core_tests.cpp)
target_link_libraries(core_tests PRIVATE epsilon_active_core)
add_test(NAME core_tests COMMAND core_tests)
