## Add libdivsufsort
set(divsufsort_cmake_file ${CMAKE_HOME_DIRECTORY}/external/libdivsufsort/CMakeLists.txt)
if(NOT EXISTS ${divsufsort_cmake_file})
    EXECUTE_PROCESS(COMMAND git submodule update --init
                    WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}
                    OUTPUT_QUIET
            )
endif(NOT EXISTS ${divsufsort_cmake_file})

set(BUILD_SHARED_LIBS OFF CACHE BOOL "Do not build a shared library for libdivsufsort")
set(BUILD_EXAMPLES OFF CACHE BOOL "Do not build libdivsufsort example")
set(BUILD_DIVSUFSORT64 ON CACHE BOOL "Build libdivsufsort in 64-bits mode")

message("CMAKE_INSTALL_INCUDEDIR=${CMAKE_INSTALL_INCUDEDIR}")    

add_subdirectory(libdivsufsort)
