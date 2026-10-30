# CMake version check
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(duneanasel VERSION 0.0.1 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Changes default install path to be a subdirectory of the build dir.
# Can set build dir at configure time with -DCMAKE_INSTALL_PREFIX=/install/path
if(CMAKE_INSTALL_PREFIX STREQUAL "" OR CMAKE_INSTALL_PREFIX STREQUAL
  "/usr/local")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")
elseif(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")
endif()

# Use the compilers found in the path
find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

option(duneanasel_build_examples "Whether to build examples" OFF)

find_package(ROOT 6 REQUIRED)

if(NOT DEFINED DUNE_ANAOBJ_BRANCH)
  SET(DUNE_ANAOBJ_BRANCH v04_02_01)
endif()

find_package(duneanaobj 4 QUIET)
if(NOT duneanaobj_FOUND)
  include(FetchContent)
  FetchContent_Declare(
      SRProxy
      GIT_REPOSITORY https://github.com/luketpickering/SRProxy.git
      GIT_TAG        feature/Clingification_rebase3
    )
  FetchContent_Declare(
      duneanaobj
      GIT_REPOSITORY https://github.com/DUNE/duneanaobj.git
      GIT_TAG        ${DUNE_ANAOBJ_BRANCH}
    )
  SET(SKIP_CET ON)
  FetchContent_MakeAvailable(SRProxy duneanaobj)
endif()

if(NOT TARGET duneanaobj::StandardRecord)
  message(FATAL_ERROR "duneanasel Expected dependency target: duneanaobj::StandardRecord")
endif()

add_library(duneanasel_all INTERFACE)

target_link_libraries(duneanasel_all INTERFACE duneanaobj::StandardRecord)
target_include_directories(duneanasel_all INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:include>)

add_library(duneanasel::all ALIAS duneanasel_all)

if(duneanasel_build_examples)
  add_subdirectory(examples)
endif()
