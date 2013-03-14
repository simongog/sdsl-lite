include_directories(#"${CMAKE_CURRENT_SOURCE_DIR}/../include"
					"${CMAKE_CURRENT_BINARY_DIR}/../include" 
					"${CMAKE_CURRENT_BINARY_DIR}/../external/libdivsufsort-2.0.1/include" 
					)


file(GLOB libFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}  "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp") # select all .cpp-files

set( sdsl_SRCS ${libFiles} )

add_library( sdsl ${sdsl_SRCS} )

install(TARGETS sdsl
		RUNTIME DESTINATION bin
		LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib)	

set_target_properties(sdsl PROPERTIES
	VERSION "${LIBRARY_VERSION_FULL}"	
	SOVERSION "${LIBRARY_VERSION_MAJOR}"	
#	DEFINE_SYMBOL SDSL_BUILD_DLL
#	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/../examples"
	)
