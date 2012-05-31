include(CheckIncludeFile)	
include(CheckIncludeFileCXX)	
include(CheckTypeSize)	

## Check for header files ##
check_include_file_CXX(cstdio HAVE_STDIO)
check_include_file(stdint.h HAVE_STDINT_H)
check_include_file(inttypes.h HAVE_INTTYPES_H)

if( HAVE_STDINT_H )	
	set( INTINCFILE "#include <stdint.h>" )
elseif( HAVE_INTTYPES_H )
	set( INTINCFILE "#include <inttypes.h>" )
else()
	set( INTINCFILE "" )
endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/uintx_t.hpp.cmake" 
		       "${CMAKE_CURRENT_BINARY_DIR}/uintx_t.hpp" @ONLY)	
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/uintx_t.hpp" DESTINATION include/sdsl)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/algorithms_for_suffix_array_construction.hpp.cmake"
		       "${CMAKE_CURRENT_BINARY_DIR}/algorithms_for_suffix_array_construction.hpp")	
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/algorithms_for_suffix_array_construction.hpp" DESTINATION include/sdsl)

## copy hpp files to the binary tree ##

file(GLOB hppFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp") # select all .hpp-files

foreach(hppFile ${hppFiles}) # copy each file
	configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/${hppFile}" "${CMAKE_CURRENT_BINARY_DIR}/${hppFile}" COPYONLY ) 
	install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${hppFile}" DESTINATION include/sdsl)
#	MESSAGE(${hppFile})
endforeach(hppFile)

