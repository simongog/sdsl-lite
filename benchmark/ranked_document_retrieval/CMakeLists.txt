cmake_minimum_required(VERSION 2.8.7)

project(ranked_document_retrieval CXX C)

INCLUDE_DIRECTORIES($ENV{HOME}/include ${CMAKE_HOME_DIRECTORY}/src)
LINK_DIRECTORIES($ENV{HOME}/lib)

set(CMAKE_CXX_FLAGS "-msse4.2 -O3 -std=c++11 -Wall -DNDEBUG")

# read the index configs
file(GLOB index_config_files RELATIVE ${CMAKE_HOME_DIRECTORY}/config/ "${CMAKE_HOME_DIRECTORY}/config/*.config")
foreach(f ${index_config_files})
  file(STRINGS ${CMAKE_HOME_DIRECTORY}/config/${f} config_contents)
  foreach(keyvalue ${config_contents})
    string(REGEX REPLACE "^[ ]+" "" keyvalue ${keyvalue})
    string(REGEX MATCH "^[^=]+" key ${keyvalue})
    string(REPLACE "${key}=" "" value ${keyvalue})
    set(${key} "${value}")
  endforeach(keyvalue)

  ADD_EXECUTABLE(build-idx-${NAME}.x src/build_idx.cpp)
  TARGET_LINK_LIBRARIES(build-idx-${NAME}.x sdsl divsufsort divsufsort64 pthread)
  set_property(TARGET build-idx-${NAME}.x PROPERTY COMPILE_DEFINITIONS IDX_TYPE="${TYPE}")

  ADD_EXECUTABLE(query-idx-${NAME}.x src/query_idx.cpp)
  TARGET_LINK_LIBRARIES(query-idx-${NAME}.x sdsl divsufsort divsufsort64 pthread)
  set_property(TARGET query-idx-${NAME}.x PROPERTY COMPILE_DEFINITIONS NAME=${NAME} IDX_TYPE="${TYPE}")

endforeach(f)