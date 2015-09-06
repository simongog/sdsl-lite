## Add libdivsufsort
set(divsufsort_cmake_file ${CMAKE_HOME_DIRECTORY}/external/libdivsufsort/CMakeLists.txt)
if(NOT EXISTS ${divsufsort_cmake_file})
    EXECUTE_PROCESS(COMMAND git submodule update --init
                    WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}
                    OUTPUT_QUIET
            )
endif(NOT EXISTS ${divsufsort_cmake_file})

set(BUILD_DIVSUFSORT64 ON)
add_subdirectory(libdivsufsort)
