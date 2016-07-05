function(find_or_download_package package package_found_prefix target_name)
    find_package(${package} ${ARGN})

    if(${package_found_prefix}_FOUND)
        message(STATUS "${package} found, no need to build locally")
    else()
        message(STATUS "${package} not found, will download and build locally")
        include(Download${package})
    endif()

    add_library(${target_name} INTERFACE)

    foreach(sublib ${${package_found_prefix}_LIBRARIES})
        get_filename_component(LIB_NAME ${sublib} NAME_WE)
        set(sublib_target ${package}_cmake_dep_${LIB_NAME})

        message(STATUS ${sublib} "  " ${sublib_target})

        # Detect wether this is a shared or a static library
        # NB: This might not work on windows
        get_filename_component(LIB_EXT ${sublib} EXT)
        if(${LIB_EXT} STREQUAL ${CMAKE_STATIC_LIBRARY_SUFFIX})
            add_library(${sublib_target} IMPORTED STATIC GLOBAL)
        else()
            add_library(${sublib_target} IMPORTED SHARED GLOBAL)
        endif()

        foreach(sublib_target_dep ${${package_found_prefix}_CMAKE_DEP})
            message(STATUS "sublib target dep: " ${sublib_target_dep})
            add_dependencies(${sublib_target} ${sublib_target_dep})
        endforeach()

        set_target_properties(${sublib_target} PROPERTIES
            "IMPORTED_LOCATION" "${sublib}"
            "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
            "INTERFACE_INCLUDE_DIRECTORIES" "${${package_found_prefix}_INCLUDE_DIRS}"
        )
        target_link_libraries(${target_name} INTERFACE
            ${sublib_target})
    endforeach(sublib)

endfunction()
