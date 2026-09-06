cmake_minimum_required(VERSION 3.20)
project(nonsmooth_gp LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
find_package(OpenMP REQUIRED)

set(LEGACY_SOURCES
  framework/legacy/src/batch_acceptance.cpp framework/legacy/src/bisection.cpp
  framework/legacy/src/bookshelf.cpp framework/legacy/src/compact_recovery.cpp
  framework/legacy/src/coarse_flow.cpp framework/legacy/src/density_coordinate.cpp
  framework/legacy/src/density.cpp framework/legacy/src/hpwl.cpp
  framework/legacy/src/lambda_controller.cpp framework/legacy/src/optimizer.cpp
  framework/legacy/src/placer.cpp framework/legacy/src/recovery.cpp
  framework/legacy/src/swap_recovery.cpp framework/legacy/src/transport.cpp)
add_library(nsgp_exact_core ${LEGACY_SOURCES})
target_include_directories(nsgp_exact_core PUBLIC framework/legacy/include framework/code third_party)
target_link_libraries(nsgp_exact_core PUBLIC OpenMP::OpenMP_CXX)
if(MSVC)
  target_compile_options(nsgp_exact_core PRIVATE /W4 /O2)
else()
  target_compile_options(nsgp_exact_core PRIVATE -Wall -Wextra -Wpedantic -O3)
endif()
add_executable(nsgp framework/code/main.cpp)
target_link_libraries(nsgp PRIVATE nsgp_exact_core)

# Historical replay tooling is intentionally separate from the default
# challenge executable. It is only used to reproduce H219/H221 in the H375
# provenance chain and must be selected explicitly by an experiment script.
add_executable(nsgp_legacy_stage framework/legacy/src/legacy_stage_main.cpp)
target_link_libraries(nsgp_legacy_stage PRIVATE nsgp_exact_core)
add_executable(nsgp_historical_homotopy
  modules/historical_dct_poisson/code/homotopy_main.cpp
  modules/historical_dct_poisson/code/electric.cpp
  modules/historical_dct_poisson/code/bookshelf.cpp
  modules/historical_dct_poisson/code/spectral.cpp)
target_include_directories(nsgp_historical_homotopy PRIVATE
  framework/legacy/include modules/historical_dct_poisson/code)
target_link_libraries(nsgp_historical_homotopy PRIVATE nsgp_exact_core OpenMP::OpenMP_CXX)
enable_testing()
add_executable(nsgp_tests tests/test_legacy_numeric_contracts.cpp)
target_link_libraries(nsgp_tests PRIVATE nsgp_exact_core)
add_test(NAME nsgp_numeric_contracts COMMAND nsgp_tests)
