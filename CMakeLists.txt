cmake_minimum_required(VERSION 2.6.0)	
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")	
include(AppendCompilerFlags)
include(CheckSSE4_2)	
	
## Project information ##
PROJECT(sdsl CXX C)
set(PROJECT_VEDOR "Simon Gog")
set(PROJECT_CONTACT "simon.gog@gmail.com")
set(PROJECT_URL "https://github.com/simongog/sdsl")
set(PROJECT_DESCRIPTION "The succinct data structure library")

## Project options ##
#option(BUILD_EXAMPLES "Build examples" OFF)
option(USE_LIBDIVSUFSORT "Use the libdivsufsort library of Yuta Mori to construct suffix arrays" ON)
if(USE_LIBDIVSUFSORT)
	find_package(divsufsort)
	find_package(divsufsort64)
else()	
	set( divsufsort_FOUND 0 )
	set( divsufsort64_FOUND 0 )
endif()	

if(NOT divsufsort_FOUND)
	set( divsufsort_INCLUDE_DIRS "." )
	message("WARNING: libdivsuftsort is not installed! ")
	message("         The build in default algorithm   ")
	message("         will use much more time and space")
	message("         then libdivsufsort!")
endif()

if(NOT divsufsort64_FOUND)
	set( divsufsort64_INCLUDE_DIRS "." )
	message("WARNING: libdivsuftsort64 is not installed! ")
	message("         The build in default algorithm   ")
	message("         will use much more time and space")
	message("         then libdivsufsort64!")
endif()	



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

#append_cxx_compiler_flags("-W -ffast-math -Wall -O9 -funroll-loops -DNDEBUG -DWRITE_R_OUTPUT -DMEM_INFO" "GCC" CMAKE_CXX_FLAGS)
append_cxx_compiler_flags("-ffast-math -Wall -pedantic -enable-long-long -Weffc++ -Werror -O9 -funroll-loops -DNDEBUG -DWRITE_R_OUTPUT -DMEM_INFO" "GCC" CMAKE_CXX_FLAGS)


if( BUILTIN_POPCNT )
	append_cxx_compiler_flags("-msse4.2" "GCC" CMAKE_CXX_FLAGS)
	message("CPU seems to support fast popcount.\n")
	message("Flag -msse4.2 added to compile options.\n")
endif()
#append_cxx_compiler_flags("-Wall -ffast-math -Wall -g -DNDEBUG" "GCC" CMAKE_CXX_FLAGS)

add_subdirectory(include)
add_subdirectory(lib)

MESSAGE("Options")
MESSAGE("USE_LIBDIVSUFSORT = ${USE_LIBDIVSUFSORT}")

## CPack configuration ##
#set(CPACK_GENERATOR "TGZ;TBZ2;ZIP")
#set(CPACK_SOURCE_GENERATOR "TGZ;TBZ2;ZIP")
#include(ProjectCPack)

