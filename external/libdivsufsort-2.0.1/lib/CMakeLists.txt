include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../include"
                    "${CMAKE_CURRENT_BINARY_DIR}/../include")

set(divsufsort_SRCS divsufsort.c sssort.c trsort.c utils.c)

## libdivsufsort ##
add_library(divsufsort ${divsufsort_SRCS})
install(TARGETS divsufsort
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
set_target_properties(divsufsort PROPERTIES
  VERSION   "${LIBRARY_VERSION_FULL}"
  SOVERSION "${LIBRARY_VERSION_MAJOR}"
  DEFINE_SYMBOL DIVSUFSORT_BUILD_DLL
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/../examples")

## libdivsufsort64 ##
if(BUILD_DIVSUFSORT64)
  add_library(divsufsort64 ${divsufsort_SRCS})
  install(TARGETS divsufsort64
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
  set_target_properties(divsufsort64 PROPERTIES
    VERSION   "${LIBRARY_VERSION_FULL}"
    SOVERSION "${LIBRARY_VERSION_MAJOR}"
    DEFINE_SYMBOL DIVSUFSORT_BUILD_DLL
    COMPILE_FLAGS "-DBUILD_DIVSUFSORT64"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/../examples")
endif(BUILD_DIVSUFSORT64)
