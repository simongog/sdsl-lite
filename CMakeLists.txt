cmake_minimum_required(VERSION 2.6.0)	
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")	
include(AppendCompilerFlags)
include(CheckSSE4_2)	
	
## Project information ##
project(sdsl CXX C)
set(PROJECT_VEDOR "Simon Gog")
set(PROJECT_CONTACT "simon.gog@gmail.com")
set(PROJECT_URL "https://github.com/simongog/sdsl-lite")
set(PROJECT_DESCRIPTION "SDSL: Succinct Data Structure Library")
set(CMAKE_BUILD_TYPE "Release")

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

append_cxx_compiler_flags("-std=c++11 -Wall -DNDEBUG" "GCC" CMAKE_CXX_FLAGS)
append_cxx_compiler_flags("-O3 -ffast-math -funroll-loops" "GCC" CMAKE_CXX_OPT_FLAGS)
append_cxx_compiler_flags("-std=c++11 -funroll-loops -DNDEBUG -stdlib=libc++" "CLANG" CMAKE_CXX_FLAGS)
append_cxx_compiler_flags("-O3 -ffast-math -funroll-loops" "CLANG" CMAKE_CXX_OPT_FLAGS)

if( BUILTIN_POPCNT )
	append_cxx_compiler_flags("-msse4.2" "GCC" CMAKE_CXX_OPT_FLAGS)
	append_cxx_compiler_flags("-msse4.2" "CLANG" CMAKE_CXX_OPT_FLAGS)
	message("CPU seems to support fast popcount.\n")
	message("Flag -msse4.2 added to compile options.\n")
endif()

add_subdirectory(external)
add_subdirectory(include)
add_subdirectory(lib)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Make.helper.cmake" 
		       "${CMAKE_CURRENT_SOURCE_DIR}/Make.helper" @ONLY)	

## Add 'uninstall' target ##
CONFIGURE_FILE(
  "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/CMakeModules/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)
ADD_CUSTOM_TARGET(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/CMakeModules/cmake_uninstall.cmake")
