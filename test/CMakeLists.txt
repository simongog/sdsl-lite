add_subdirectory(gtest-1.6.0)
enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

include_directories("${CMAKE_CURRENT_BINARY_DIR}/../include"   # for sdsl 
					${divsufsort_INCLUDE_DIRS}               # for divsufsort
					${divsufsort64_INCLUDE_DIRS}             # for divsufsort64
					"${gtest_SOURCE_DIR}/include"
					) 

link_directories("${CMAKE_INSTALL_PREFIX}/lib")

file(GLOB test_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/*Test.cpp") # select all .cpp-files

foreach(test_SRC ${test_SRCS})
	string(REPLACE ".cpp" "" test_PROG ${test_SRC})
	message("...${test_PROG} ${test_SRC}")
	add_executable(${test_PROG} ${test_SRC})
	target_link_libraries(${test_PROG} gtest)
	target_link_libraries(${test_PROG} sdsl)
	target_link_libraries(${test_PROG} divsufsort)
	target_link_libraries(${test_PROG} divsufsort64)
	add_test(${test_PROG} ${test_PROG})
endforeach(testSource)

