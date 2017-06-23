## Add libdivsufsort
set(divsufsort_cmake_file ${CMAKE_CURRENT_LIST_DIR}/libdivsufsort/CMakeLists.txt)
if(NOT EXISTS ${divsufsort_cmake_file})
    EXECUTE_PROCESS(COMMAND git submodule update --init external/libdivsufsort
                    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/..
                    OUTPUT_QUIET
            )
endif(NOT EXISTS ${divsufsort_cmake_file})

set(BUILD_SHARED_LIBS OFF CACHE BOOL "Do not build a shared library for libdivsufsort")
set(BUILD_EXAMPLES OFF CACHE BOOL "Do not build libdivsufsort example")
set(BUILD_DIVSUFSORT64 ON CACHE BOOL "Build libdivsufsort in 64-bits mode")

add_subdirectory(libdivsufsort)

## Add gtest
set(gtest_dir ${CMAKE_CURRENT_LIST_DIR}/googletest)
set(gtest_file ${gtest_dir}/CMakeLists.txt)
if(NOT EXISTS ${gtest_file})
    EXECUTE_PROCESS(COMMAND git submodule update --init external/googletest
                    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/..
                    OUTPUT_QUIET
            )
endif(NOT EXISTS ${gtest_file})

add_subdirectory(googletest/googletest)


