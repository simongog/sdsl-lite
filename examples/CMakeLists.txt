file(GLOB cppFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp") # select all .hpp-files


set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

include(ExternalProject)
include(find_or_download_package)

message("Boost root set to: ${BOOST_ROOT}")

find_or_download_package(Boost Boost boost COMPONENTS system filesystem regex)

INCLUDE_DIRECTORIES(
        ${CMAKE_CURRENT_BINARY_DIR}/../external/libdivsufsort/include
        ${Boost_INCLUDE_DIRS}
        ${gtest_SOURCE_DIR}/include
        ${STXXL_INCLUDE_DIRS}
        )

LINK_LIBRARIES(sdsl gtest divsufsort divsufsort64 stxxl boost)

# Set compiler flags to show all warnings (MSVC and GCC)
if(MSVC)
    # Force to always compile with W4
    if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
        string(REGEX REPLACE "/W[0-4]" "/W3" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3")
    endif()
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    # Update if necessary
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wall -std=c++11 -fsanitize=memory -ltbb -funroll-loops -Wno-long-long -pedantic -Wno-write-strings -flto=10 -g -pipe" )

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -fopenmp -Wall -std=c++11 -Wno-long-long -pedantic -Wno-write-strings -g -pipe -D_GLIBCXX_USE_CXX11_ABI=0" )
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -funroll-loops -O3 -flto=10 ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -D_GLIBCXX_DEBUG -DDEBUG -ggdb -fsanitize=address") #great options for weird segfaults  -fsanitize=undefined
    # disable this by default
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -flto=10" )
endif()

foreach(cppFile ${cppFiles}) # copy each file
    GET_FILENAME_COMPONENT(cpp_name_we ${cppFile} NAME_WE)
    ADD_EXECUTABLE(${cpp_name_we} ${cpp_name_we}.cpp)
endforeach(cppFile)
