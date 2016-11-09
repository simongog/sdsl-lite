file(GLOB hFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/dacs/*.h") # select all .h-files

foreach(hFile ${hFiles}) # copy each file
    configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/${hFile}" "${CMAKE_CURRENT_BINARY_DIR}/${hFile}" COPYONLY )
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${hFile}" DESTINATION include/dacs)
    #	MESSAGE(${hFile})
endforeach(hFile)
