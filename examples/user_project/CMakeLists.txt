# CMake version check
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(duneanauserexample VERSION 0.0.1 LANGUAGES CXX)

 include(FetchContent)

FetchContent_Declare(
  duneanasel
  GIT_REPOSITORY https://github.com/luketpickering/duneanasel.git
)

FetchContent_MakeAvailable(duneanasel)

add_executable(myana myana.cxx)
target_link_libraries(myana duneanasel::all)
