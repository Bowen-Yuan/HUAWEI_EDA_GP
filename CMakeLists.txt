cmake_minimum_required(VERSION 3.20)
project(nonsmooth_gp LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
find_package(OpenMP REQUIRED)

set(NUMERIC_KERNEL_SOURCES
  framework/kernel/src/bookshelf.cpp
  framework/kernel/src/density.cpp
  framework/kernel/src/hpwl.cpp
  framework/kernel/src/optimizer.cpp
  framework/kernel/src/step_policy.cpp
  framework/kernel/src/lambda_controller.cpp
  framework/kernel/src/nonlocal_descent.cpp)
set(MODULE_ALGORITHM_SOURCES
  modules/exact_joint_gp/code/batch_acceptance.cpp
  modules/exact_joint_gp/code/placer.cpp
  modules/exact_recovery/code/compact_recovery.cpp
  modules/exact_recovery/code/recovery.cpp
  modules/surplus_bisection/code/bisection.cpp
  modules/equal_shape_swap/code/swap_recovery.cpp
  modules/global_capacity_transport/code/coarse_flow.cpp
  modules/global_capacity_transport/code/transport.cpp
  modules/density_coordinate/code/density_coordinate.cpp)
add_library(nsgp_numeric_kernel ${NUMERIC_KERNEL_SOURCES} ${MODULE_ALGORITHM_SOURCES})
target_include_directories(nsgp_numeric_kernel PUBLIC framework/kernel/include framework/code third_party)
target_link_libraries(nsgp_numeric_kernel PUBLIC OpenMP::OpenMP_CXX)
if(MSVC)
  target_compile_options(nsgp_numeric_kernel PRIVATE /W4 /O2)
else()
  target_compile_options(nsgp_numeric_kernel PRIVATE -Wall -Wextra -Wpedantic -O3)
endif()
add_executable(nsgp
  framework/code/main.cpp
  framework/code/microkernel.cpp
  modules/global_view_gp/code/global_view_lab.cpp
  modules/layout_init/code/module.cpp
  modules/hpwl_adam/code/module.cpp
  modules/exact_joint_gp/code/module.cpp
  modules/exact_recovery/code/module.cpp
  modules/surplus_bisection/code/module.cpp
  modules/equal_shape_swap/code/module.cpp)
target_sources(nsgp PRIVATE
  modules/global_capacity_transport/code/module.cpp
  modules/density_coordinate/code/module.cpp)
target_include_directories(nsgp PRIVATE modules/global_view_gp/code)
target_sources(nsgp PRIVATE
  modules/adaptive_lambda_gp/code/adaptive_lambda_gp.cpp
  modules/adaptive_lambda_gp/code/module.cpp)
target_sources(nsgp PRIVATE
  modules/density_multiscale_active_gp/code/module.cpp
  modules/density_cut_pressure_gp/code/module.cpp
  modules/density_transport_gp/code/module.cpp
  modules/density_charge_gp/code/module.cpp
  modules/finite_radius_oracle_gp/code/module.cpp)
target_include_directories(nsgp PRIVATE modules/adaptive_lambda_gp/code)
target_include_directories(nsgp PRIVATE modules)
target_link_libraries(nsgp PRIVATE nsgp_numeric_kernel)
if(WIN32)
  target_link_libraries(nsgp PRIVATE Advapi32)
endif()
execute_process(COMMAND git rev-parse --short HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE NSGP_GIT_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)
execute_process(COMMAND git branch --show-current
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE NSGP_GIT_BRANCH OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)
execute_process(COMMAND git remote get-url origin
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE NSGP_GIT_REMOTE OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)
execute_process(COMMAND git status --porcelain --untracked-files=no
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE NSGP_GIT_STATUS OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)
if(NSGP_GIT_STATUS STREQUAL "")
  set(NSGP_GIT_DIRTY "clean")
else()
  set(NSGP_GIT_DIRTY "dirty")
endif()
target_compile_definitions(nsgp PRIVATE
  NSGP_GIT_COMMIT="${NSGP_GIT_COMMIT}"
  NSGP_GIT_BRANCH="${NSGP_GIT_BRANCH}"
  NSGP_GIT_DIRTY="${NSGP_GIT_DIRTY}"
  NSGP_GIT_REMOTE="${NSGP_GIT_REMOTE}")

# Historical replay tooling is intentionally separate from the default
# challenge executable. It is only used to reproduce H219/H221 in the H375
# provenance chain and must be selected explicitly by an experiment script.
add_executable(nsgp_legacy_stage modules/historical_exact_replay/code/legacy_stage_main.cpp)
target_link_libraries(nsgp_legacy_stage PRIVATE nsgp_numeric_kernel)
add_executable(nsgp_historical_homotopy
  modules/historical_dct_poisson/code/homotopy_main.cpp
  modules/historical_dct_poisson/code/electric.cpp
  modules/historical_dct_poisson/code/bookshelf.cpp
  modules/historical_dct_poisson/code/spectral.cpp)
target_include_directories(nsgp_historical_homotopy PRIVATE
  framework/kernel/include modules/historical_dct_poisson/code)
target_link_libraries(nsgp_historical_homotopy PRIVATE nsgp_numeric_kernel OpenMP::OpenMP_CXX)
enable_testing()
add_executable(nsgp_tests tests/test_numeric_contracts.cpp)
target_link_libraries(nsgp_tests PRIVATE nsgp_numeric_kernel)
add_test(NAME nsgp_numeric_contracts COMMAND nsgp_tests)

add_executable(nsgp_adaptive_tests tests/test_adaptive_lambda_contracts.cpp)
target_sources(nsgp_adaptive_tests PRIVATE
  modules/adaptive_lambda_gp/code/adaptive_lambda_gp.cpp)
target_link_libraries(nsgp_adaptive_tests PRIVATE nsgp_numeric_kernel)
target_include_directories(nsgp_adaptive_tests PRIVATE modules/adaptive_lambda_gp/code)
add_test(NAME nsgp_adaptive_lambda_contracts COMMAND nsgp_adaptive_tests)

add_executable(nsgp_nonlocal_oracle_tests tests/test_nonlocal_oracle_contracts.cpp)
target_link_libraries(nsgp_nonlocal_oracle_tests PRIVATE nsgp_numeric_kernel)
add_test(NAME nsgp_nonlocal_oracle_contracts COMMAND nsgp_nonlocal_oracle_tests)
