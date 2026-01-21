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

# download CPM.cmake
file(
  DOWNLOAD
  https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.40.2/CPM.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
)
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)

option(duneanasel_USE_SRProxy "Whether to build proxy classes for Standard Record objects" OFF)
option(duneanasel_SRProxy_PYTHON_ENABLED "Whether to build python bindings for proxy classes" OFF)

find_package(ROOT 6 REQUIRED)

find_package(duneanaobj QUIET)

if(duneanaobj_FOUND AND TARGET duneanaobj::StandardRecordProxy)
  message(STATUS "Found duneanaobj in ${duneanaobj_SOURCE}.")
  add_library(duneanaobj_all INTERFACE)
  target_link_libraries(duneanaobj_all PUBLIC duneanaobj::StandardRecordProxy duneanaobj::StandardRecordFlat)
  add_library(duneanaobj::all ALIAS duneanaobj_all)
elseif(duneanasel_USE_SRProxy)
  if(NOT duneanaobj_FOUND)
    CPMFindPackage(
      NAME duneanaobj
      GIT_TAG feature/Clingification_rebase2024
      GITHUB_REPOSITORY luketpickering/duneanaobj
      VERSION 4.0.0
      OPTIONS
        "duneanaobj_PYTHON_ENABLED" "${duneanasel_SRProxy_PYTHON_ENABLED}"
    )
  endif()

  if(NOT TARGET duneanaobj::all)
    message(FATAL_ERROR "duneanasel Expected dependency target: duneanaobj::all")
  endif()
else()
  if(NOT DEFINED DUNE_ANAOBJ_BRANCH OR "${DUNE_ANAOBJ_BRANCH}x" STREQUAL "x")
    set(DUNE_ANAOBJ_BRANCH v03_06_01)
  endif()
  CPMAddPackage(
    NAME duneanaobj
    GIT_TAG ${DUNE_ANAOBJ_BRANCH}
    GITHUB_REPOSITORY DUNE/duneanaobj
    DOWNLOAD_ONLY YES
  )
  include_directories(${duneanaobj_SOURCE_DIR})

  ROOT_GENERATE_DICTIONARY(StandardRecordDict
    ${duneanaobj_SOURCE_DIR}/duneanaobj/StandardRecord/StandardRecord.h
    LINKDEF
      ${duneanaobj_SOURCE_DIR}/duneanaobj/StandardRecord/classes_def.xml)

  file(GLOB SR_IMPL_FILES "${duneanaobj_SOURCE_DIR}/duneanaobj/StandardRecord/*.cxx")
  LIST(APPEND SR_IMPL_FILES ${CMAKE_CURRENT_BINARY_DIR}/StandardRecordDict.cxx)

  add_library(duneanaobj_StandardRecord SHARED ${SR_IMPL_FILES})
  target_link_libraries(duneanaobj_StandardRecord PUBLIC ROOT::MathCore)

  target_include_directories(duneanaobj_StandardRecord PUBLIC
    $<BUILD_INTERFACE:${duneanaobj_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>)
  target_include_directories(duneanaobj_StandardRecord PRIVATE
    $<BUILD_INTERFACE:${duneanaobj_SOURCE_DIR}/duneanaobj/StandardRecord> #root puts this in the dictionary
    )

  set_target_properties(duneanaobj_StandardRecord PROPERTIES EXPORT_NAME all)
  install(TARGETS duneanaobj_StandardRecord EXPORT duneanasel-targets DESTINATION lib)

  file(GLOB SR_HEADER_FILES "${duneanaobj_SOURCE_DIR}/duneanaobj/StandardRecord/*.h")

  install(FILES ${SR_HEADER_FILES} DESTINATION include/duneanaobj/StandardRecord)
  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/libStandardRecordDict_rdict.pcm
    ${CMAKE_CURRENT_BINARY_DIR}/libStandardRecordDict.rootmap
    DESTINATION lib)

  add_library(duneanaobj::all ALIAS duneanaobj_StandardRecord)

endif()

add_library(duneanasel_all INTERFACE)

target_link_libraries(duneanasel_all INTERFACE duneanaobj::all)
target_include_directories(duneanasel_all INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
  $<INSTALL_INTERFACE:include>)

add_library(duneanasel::all ALIAS duneanasel_all)

add_subdirectory(examples)
