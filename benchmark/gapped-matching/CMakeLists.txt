cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0015 NEW)
set(CMAKE_MODULE_PATH "${CMAKE_HOME_DIRECTORY}/../../CMakeModules")
include(AppendCompilerFlags)
include(ExternalProject)

project(gapped-matching CXX C)
set(CMAKE_BUILD_TYPE "Release")

# parse make helper (HACK  :( )
file(STRINGS ${CMAKE_HOME_DIRECTORY}/../../Make.helper makehelper_contents)
list (GET makehelper_contents 0 INC_DIR_STR)
list (GET makehelper_contents 1 LIB_DIR_STR)
STRING(REGEX REPLACE "INC_DIR = (.*)" "\\1" SDSL_INC_DIR "${LIB_DIR_STR}")
STRING(REGEX REPLACE "LIB_DIR = (.*)" "\\1" SDSL_LIB_DIR "${INC_DIR_STR}")
MESSAGE( STATUS "SDSL_INC_DIR: " ${SDSL_INC_DIR} )
MESSAGE( STATUS "SDSL_LIB_DIR: " ${SDSL_LIB_DIR} )


# set include and lib dirs
INCLUDE_DIRECTORIES(${SDSL_INC_DIR} 
					${CMAKE_HOME_DIRECTORY}/external/easyloggingpp/src/
					${CMAKE_HOME_DIRECTORY}/include
					)
LINK_DIRECTORIES(${SDSL_LIB_DIR})


# C++11 compiler Check
if(NOT CMAKE_CXX_COMPILER_VERSION) # work around for cmake versions smaller than 2.8.10
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE CMAKE_CXX_COMPILER_VERSION)
endif()
if(CMAKE_CXX_COMPILER MATCHES ".*clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      set(CMAKE_COMPILER_IS_CLANGXX 1)
endif()
if( (CMAKE_COMPILER_IS_GNUCXX AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 4.7) OR
    (CMAKE_COMPILER_IS_CLANGXX AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 3.2))
  message(FATAL_ERROR "Your C++ compiler does not support C++11. Please install g++ 4.7 (or greater) or clang 3.2 (or greater)")
else()
  message(STATUS "Compiler is recent enough to support C++11.")
endif()
if( MEASURE_ENERGY )
    message(STATUS "Measure energy.")
    append_cxx_compiler_flags("-DLIKWID_PERFMON -L$ENV{LIKWID_LIB} -I$ENV{LIKWID_INCLUDE} -llikwid" "GCC" CMAKE_CXX_FLAGS)
endif()
if( CMAKE_COMPILER_IS_GNUCXX )
    append_cxx_compiler_flags("-std=c++11 -Wall -Wextra " "GCC" CMAKE_CXX_FLAGS)
    append_cxx_compiler_flags("-msse4.2 -O3 -ffast-math -funroll-loops" "GCC" CMAKE_CXX_FLAGS_RELEASE)
else()
    append_cxx_compiler_flags("-std=c++11" "CLANG" CMAKE_CXX_FLAGS)
    append_cxx_compiler_flags("-stdlib=libc++" "CLANG" CMAKE_CXX_FLAGS)
    append_cxx_compiler_flags("-msse4.2 -O3  -ffast-math -funroll-loops -DNDEBUG" "CLANG" CMAKE_CXX_FLAGS_RELEASE)
endif()

FIND_PACKAGE(Boost 1.45 COMPONENTS regex)

# # read the index configs
file(GLOB index_config_files RELATIVE ${CMAKE_HOME_DIRECTORY}/config/ "${CMAKE_HOME_DIRECTORY}/config/*.config")
foreach(f ${index_config_files})
	file(STRINGS ${CMAKE_HOME_DIRECTORY}/config/${f} config_contents)
	set(compile_defs "")
	foreach(keyvalue ${config_contents})
		string(REGEX REPLACE "^[ ]+" "" keyvalue ${keyvalue})
		string(REGEX MATCH "^[^=]+" key ${keyvalue})
		string(REPLACE "${key}=" "" value ${keyvalue})
		set(${key} "${value}")
		list(APPEND compile_defs ${key}=${value})
	endforeach(keyvalue)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
	ADD_EXECUTABLE(gm_index-${NAME}.x src/gm_index.cpp)
	TARGET_LINK_LIBRARIES(gm_index-${NAME}.x sdsl divsufsort divsufsort64 pthread ${Boost_LIBRARIES})
	set_property(TARGET gm_index-${NAME}.x PROPERTY COMPILE_DEFINITIONS IDXNAME="${NAME}" ${compile_defs})

    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
	ADD_EXECUTABLE(gm_search-${NAME}.x src/gm_search.cpp)
	TARGET_LINK_LIBRARIES(gm_search-${NAME}.x sdsl divsufsort divsufsort64 pthread ${Boost_LIBRARIES})
	set_property(TARGET gm_search-${NAME}.x PROPERTY COMPILE_DEFINITIONS IDXNAME="${NAME}" ${compile_defs})
endforeach(f)

ADD_EXECUTABLE(create_collection.x src/create_collection.cpp)
TARGET_LINK_LIBRARIES(create_collection.x sdsl divsufsort divsufsort64 pthread)

