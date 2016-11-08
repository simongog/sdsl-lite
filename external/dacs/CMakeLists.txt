file(GLOB_RECURSE DACS_SRC_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cc *.c)

add_subdirectory(include)

add_library(dacs SHARED ${DACS_SRC_FILES})

install(FILES
  "LICENSE"
  DESTINATION .
)

install(TARGETS dacs
  LIBRARY DESTINATION lib
  COMPONENT libraries
)

set(DACS_INCLUDE_DIR, ${DACS_SRC_FILES}/../include)