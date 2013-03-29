configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Make.helper.cmake" 
		       "${CMAKE_CURRENT_SOURCE_DIR}/Make.helper" @ONLY)	
add_subdirectory(indexing_count)
