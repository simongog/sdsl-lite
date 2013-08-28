
file(GLOB gtestFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
foreach(gtestFile ${gtestFiles}) 
     configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/${gtestFile}" "${CMAKE_CURRENT_BINARY_DIR}/${gtestFile}" COPYONLY )
#     install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${gtestFile}" DESTINATION include/gtest)
endforeach(gtestFile)

add_subdirectory(internal)
