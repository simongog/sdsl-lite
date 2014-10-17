
INCLUDE_DIRECTORIES(${CMAKE_HOME_DIRECTORY}/include
                    ${CMAKE_HOME_DIRECTORY}/lecture/external/cppzmq
                    ${CMAKE_HOME_DIRECTORY}/lecture/external/zeromq4-1/include
                    ${CMAKE_BINARY_DIR}/external/libdivsufsort-2.0.1/include
)

LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lecture/external/zeromq4-1/lib)

ADD_SUBDIRECTORY(external/zeromq4-1)
SET_PROPERTY(DIRECTORY external/zeromq4-1 PROPERTY ZMQ_BUILD_TESTS FALSE)

ADD_EXECUTABLE(index-daemon index-daemon.cpp)
TARGET_LINK_LIBRARIES(index-daemon sdsl divsufsort divsufsort64 libzmq)

ADD_EXECUTABLE(index-client index-client.cpp)
TARGET_LINK_LIBRARIES(index-client sdsl divsufsort divsufsort64 libzmq)
