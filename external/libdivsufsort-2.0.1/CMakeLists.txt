### cmake file for building libdivsufsort Package ###
cmake_minimum_required(VERSION 2.4)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")
include(AppendCompilerFlags)

## SVN revision ##
set(SVN_REVISION "")
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.svn")
  execute_process(COMMAND svn info --xml
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE SVN_INFO ERROR_QUIET)
  if(SVN_INFO)
    string(REGEX MATCH "<entry[^>]+" SVN_REVISION "${SVN_INFO}")
    string(REGEX REPLACE "^.*revision=\"([0-9]+)\".*$" "\\1" SVN_REVISION "${SVN_REVISION}")
  endif(SVN_INFO)
endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.svn")

## Project information ##
project(libdivsufsort C)
set(PROJECT_VENDOR "Yuta Mori")
set(PROJECT_CONTACT "yuta.256@gmail.com")
set(PROJECT_URL "http://libdivsufsort.googlecode.com/")
set(PROJECT_DESCRIPTION "A lightweight suffix sorting library")
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" PROJECT_VERSION_FULL)
string(REGEX REPLACE "[\n\r]" "" PROJECT_VERSION_FULL "${PROJECT_VERSION_FULL}")
string(REGEX REPLACE "^([0-9]+)\\.[0-9]+\\.[0-9]+$" "\\1" PROJECT_VERSION_MAJOR "${PROJECT_VERSION_FULL}")
string(REGEX REPLACE "^[0-9]+\\.([0-9]+)\\.[0-9]+$" "\\1" PROJECT_VERSION_MINOR "${PROJECT_VERSION_FULL}")
string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+)$" "\\1" PROJECT_VERSION_PATCH "${PROJECT_VERSION_FULL}")
set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}")
math(EXPR LIBRARY_VERSION_MAJOR "1 + ${PROJECT_VERSION_MAJOR}")
set(LIBRARY_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
set(LIBRARY_VERSION_PATCH "${PROJECT_VERSION_PATCH}")
set(LIBRARY_VERSION "${LIBRARY_VERSION_MAJOR}.${LIBRARY_VERSION_MINOR}")
set(LIBRARY_VERSION_FULL "${LIBRARY_VERSION}.${LIBRARY_VERSION_PATCH}")
if(SVN_REVISION)
  set(PROJECT_VERSION_FULL "svn-r${SVN_REVISION}")
endif(SVN_REVISION)

## CPack configuration ##
set(CPACK_GENERATOR "TGZ;TBZ2;ZIP")
set(CPACK_SOURCE_GENERATOR "TGZ;TBZ2;ZIP")
include(ProjectCPack)

## Project options ##
option(BUILD_SHARED_LIBS "Set to OFF to build static libraries" OFF)
option(BUILD_EXAMPLES "Build examples" OFF)
option(BUILD_DIVSUFSORT64 "Build libdivsufsort64" ON)
option(USE_OPENMP "Use OpenMP for parallelization" OFF)
option(WITH_LFS "Enable Large File Support" ON)

## Build type ##
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_VERBOSE_MAKEFILE ON)
endif(NOT CMAKE_BUILD_TYPE)

## Compiler options ##
if(MSVC)
  append_c_compiler_flags("/W4" "VC" CMAKE_C_FLAGS)
  append_c_compiler_flags("/Oi;/Ot;/Ox;/Oy" "VC" CMAKE_C_FLAGS_RELEASE)
  if(USE_OPENMP)
    append_c_compiler_flags("/openmp" "VC" CMAKE_C_FLAGS)
  endif(USE_OPENMP)
elseif(BORLAND)
  append_c_compiler_flags("-w" "BCC" CMAKE_C_FLAGS)
  append_c_compiler_flags("-Oi;-Og;-Os;-Ov;-Ox" "BCC" CMAKE_C_FLAGS_RELEASE)
else(MSVC)
  if(CMAKE_COMPILER_IS_GNUCC)
    append_c_compiler_flags("-Wall" "GCC" CMAKE_C_FLAGS)
    append_c_compiler_flags("-fomit-frame-pointer" "GCC" CMAKE_C_FLAGS_RELEASE)
    if(USE_OPENMP)
      append_c_compiler_flags("-fopenmp" "GCC" CMAKE_C_FLAGS)
    endif(USE_OPENMP)
  else(CMAKE_COMPILER_IS_GNUCC)
    append_c_compiler_flags("-Wall" "UNKNOWN" CMAKE_C_FLAGS)
    append_c_compiler_flags("-fomit-frame-pointer" "UNKNOWN" CMAKE_C_FLAGS_RELEASE)
    if(USE_OPENMP)
      append_c_compiler_flags("-fopenmp;-openmp;-omp" "UNKNOWN" CMAKE_C_FLAGS)
    endif(USE_OPENMP)
  endif(CMAKE_COMPILER_IS_GNUCC)
endif(MSVC)

## Add definitions ##
add_definitions(-DHAVE_CONFIG_H=1 -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -D__STDC_FORMAT_MACROS)

## Add subdirectories ##
#add_subdirectory(pkgconfig)
add_subdirectory(include)
add_subdirectory(lib)
if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif(BUILD_EXAMPLES)
