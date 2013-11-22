cmake_minimum_required(VERSION 2.8.7)

# compile stuff
project(ranked_document_retrieval CXX C)
include_directories($ENV{HOME}/include ${CMAKE_HOME_DIRECTORY}/src)
link_directories($ENV{HOME}/lib)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2 -O3 -std=c++11 -Wall -DNDEBUG")

# add misc binaries
add_executable(create_docweights_len.x src/create_docweights_len.cpp)
target_link_libraries(create_docweights_len.x sdsl divsufsort divsufsort64 pthread)


# read the index configs
add_custom_target(create_indexes)
add_custom_target(create_indexes_int)
file(GLOB index_config_files RELATIVE ${CMAKE_HOME_DIRECTORY}/config/ "${CMAKE_HOME_DIRECTORY}/config/*.index")
foreach(f ${index_config_files})
  file(STRINGS ${CMAKE_HOME_DIRECTORY}/config/${f} config_contents)
  foreach(keyvalue ${config_contents})
    string(REGEX REPLACE "^[ ]+" "" keyvalue ${keyvalue})
    string(REGEX MATCH "^[^=]+" key ${keyvalue})
    string(REPLACE "${key}=" "" value ${keyvalue})
    set(${key} "${value}")
  endforeach(keyvalue)

  add_executable(build-idx-${NAME}.x src/build_idx.cpp)
  target_link_libraries(build-idx-${NAME}.x sdsl divsufsort divsufsort64 pthread)
  set_property(TARGET build-idx-${NAME}.x PROPERTY COMPILE_DEFINITIONS NAME="${NAME}" IDX_TYPE=${TYPE})
  add_dependencies(create_indexes build-idx-${NAME}.x)

  add_executable(build-idx-${NAME}-int.x src/build_idx.cpp)
  target_link_libraries(build-idx-${NAME}-int.x sdsl divsufsort divsufsort64 pthread)
  set_property(TARGET build-idx-${NAME}-int.x PROPERTY COMPILE_DEFINITIONS NAME="${NAME}INT" IDX_TYPE=${TYPEINT})
  add_dependencies(create_indexes_int build-idx-${NAME}-int.x)

  add_executable(query-idx-${NAME}.x src/query_idx.cpp)
  target_link_libraries(query-idx-${NAME}.x sdsl divsufsort divsufsort64 pthread)
  set_property(TARGET query-idx-${NAME}.x PROPERTY COMPILE_DEFINITIONS NAME="${NAME}INT" IDX_TYPE=${TYPEINT})

  add_executable(query-idx-${NAME}-int.x src/query_idx.cpp)
  target_link_libraries(query-idx-${NAME}-int.x sdsl divsufsort divsufsort64 pthread)
  set_property(TARGET query-idx-${NAME}-int.x PROPERTY COMPILE_DEFINITIONS NAME="${NAME}INT" IDX_TYPE=${TYPEINT})
endforeach(f)

# read the collection configs
file(GLOB collection_config_files RELATIVE ${CMAKE_HOME_DIRECTORY}/config/ "${CMAKE_HOME_DIRECTORY}/config/*.collection")
foreach(f ${collection_config_files})
  file(STRINGS ${CMAKE_HOME_DIRECTORY}/config/${f} collection_contents)
  foreach(keyvalue ${collection_contents})
    string(REGEX REPLACE "^[ ]+" "" keyvalue ${keyvalue})
    string(REGEX MATCH "^[^=]+" key ${keyvalue})
    string(REPLACE "${key}=" "" value ${keyvalue})
    set(${key} "${value}")
  endforeach(keyvalue)

  add_custom_target(collection_${NAME} ALL DEPENDS create_indexes COMMAND ${CMAKE_HOME_DIRECTORY}/build_indexes.sh ${NAME} ${FILE} ${DOCWEIGHTS})
  add_custom_target(collection_${NAME}_int ALL DEPENDS create_indexes_int COMMAND ${CMAKE_HOME_DIRECTORY}/build_indexes_int.sh ${NAME} ${FILEINT} ${DOCWEIGHTS})

endforeach(f)
