file(GLOB cppFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp") # select all .hpp-files


set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost COMPONENTS system filesystem regex tokenizer)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../include
        ${CMAKE_CURRENT_BINARY_DIR}/../external/libdivsufsort/include
        ${Boost_INCLUDE_DIRS}
        ${gtest_SOURCE_DIR}/include
        )

LINK_LIBRARIES(sdsl gtest divsufsort divsufsort64)

foreach(cppFile ${cppFiles}) # copy each file
    GET_FILENAME_COMPONENT(cpp_name_we ${cppFile} NAME_WE)
    ADD_EXECUTABLE(${cpp_name_we} ${cpp_name_we}.cpp)
endforeach(cppFile)