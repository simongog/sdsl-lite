include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Makefile.cmake" 
		       "${CMAKE_CURRENT_SOURCE_DIR}/Makefile" @ONLY)
