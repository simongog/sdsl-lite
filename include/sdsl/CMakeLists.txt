include(CheckIncludeFile)	
include(CheckIncludeFileCXX)	
include(CheckTypeSize)	

## Check for header files ##
check_include_file_CXX(cstdio HAVE_STDIO)

## copy hpp files to the binary tree ##
include_directories(/usr/local/include/c++/4.9.3)

message(STATUS "In File examples cmake current_source_dir: ${CMAKE_CURRENT_SOURCE_DIR}")
file(GLOB hppFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp") # select all .hpp-files

foreach(hppFile ${hppFiles}) # copy each file
	configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/${hppFile}" "${CMAKE_CURRENT_BINARY_DIR}/${hppFile}" COPYONLY ) 
	install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${hppFile}" DESTINATION include/sdsl)
#	MESSAGE(${hppFile})
endforeach(hppFile)

