cmake_minimum_required(VERSION 3.8)
project(libirc)

# + ----------------- +
# | BUILDING SETTINGS |
# + ----------------- +

if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# + ----------------- +
# | LANGUAGE SETTINGS |
# + ----------------- +

# C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Flags
set(CMAKE_CXX_FLAGS -Wall)

# + ----------- +
# | DIRECTORIES |
# + ----------- +

set(IRC_MOLECULES_DIR ${CMAKE_SOURCE_DIR}/molecules/)

configure_file( ${CMAKE_SOURCE_DIR}/src/config.h.in
                ${CMAKE_BINARY_DIR}/src/config.h)

include_directories( include )
include_directories( include/libirc )
include_directories( include/libirc/wrapper )
include_directories( ${CMAKE_BINARY_DIR}/src)

# + ----------- +
# | IRC LIBRARY |
# + ----------- +

add_library( irc SHARED
                 src/atoms/atom.cpp
                 src/atoms/periodic_table.cpp
                 src/tools/mathtools.cpp
        )

set_target_properties(irc PROPERTIES LINKER_LANGUAGE CXX)

# + ----------------------- +
# | FIND EXTERNAL LIBRARIES |
# + ----------------------- +

find_package(
        Boost
        1.65.0
        REQUIRED
        )
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(irc ${Boost_LIBRARIES})

find_package(
        Armadillo
        8.200.0
        REQUIRED
        )
if( ${ARMADILLO_FOUND})
        include_directories(${ARMADILLO_INCLUDE_DIRS})
        target_link_libraries(irc ${ARMADILLO_LIBRARIES})
        add_definitions(-DHAVE_ARMA)

        message("-- Armadillo version: " ${ARMADILLO_VERSION_STRING})
endif()

#find_package(
#        Eigen3
#        3.2
#)
#if( ${EIGEN3_FOUND})
#        include_directories(${EIGEN3_INCLUDE_DIRS})
#        target_link_libraries(libirc ${EIGEN3_LIBRARIES})
#        add_definitions(-DHAVE_EIGEN3)

#        message("-- Eigen3 version: " ${EIGEN3_VERSION})
#endif()

# + ------- +
# | INSTALL |
# + ------- +

install(DIRECTORY include/libirc DESTINATION include)
install(FILES include/libirc.h DESTINATION include)
install(FILES "${CMAKE_BINARY_DIR}/src/config.h" DESTINATION include/libirc)
install(TARGETS irc DESTINATION lib)

# + ------- +
# | TESTING |
# + ------- +

enable_testing()
add_subdirectory(src/test)