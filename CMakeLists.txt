cmake_minimum_required(VERSION 3.8)
project(irc)

# + ----------------- +
# | LANGUAGE SETTINGS |
# + ----------------- +

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# + ----------- +
# | DIRECTORIES |
# + ----------- +

set(IRC_MOLECULES_DIR ${CMAKE_SOURCE_DIR}/test/)

configure_file( ${CMAKE_SOURCE_DIR}/src/config.h.in
                ${CMAKE_BINARY_DIR}/src/config.h)

include_directories( src )
include_directories( ${CMAKE_BINARY_DIR}/src)

# + ----------- +
# | IRC LIBRARY |
# + ----------- +

add_library( libirc SHARED
                    src/atoms/atom.cpp
                    src/atoms/periodic_table.cpp
                    src/atoms/periodic_table.cpp
                    src/tools/math.cpp
        )

set_target_properties(libirc PROPERTIES OUTPUT_NAME "irc")
set_target_properties(libirc PROPERTIES LINKER_LANGUAGE CXX)

# + ----------------------- +
# | FIND EXTERNAL LIBRARIES |
# + ----------------------- +

find_package(
        Boost
        1.65.0
        REQUIRED
        )
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(libirc ${Boost_LIBRARIES})

find_package(
        Armadillo
        8.200.0
        REQUIRED
        )
if( ${ARMADILLO_FOUND})
        include_directories(${ARMADILLO_INCLUDE_DIRS})
        target_link_libraries(libirc ${ARMADILLO_LIBRARIES})
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
# | TESTING |
# + ------- +

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_testing()

add_executable( test_math src/tools/math_test.cpp )
target_link_libraries( test_math libirc )
add_test( UNITTEST_math bin/test_math )

add_executable( test_periodic_table src/atoms/periodic_table_test.cpp )
target_link_libraries( test_periodic_table libirc )
add_test( UNITTEST_periodic_table bin/test_periodic_table )

add_executable( test_atom src/atoms/atom_test.cpp )
target_link_libraries( test_atom libirc )
add_test( UNITTEST_atom bin/test_atom )

add_executable( test_molecule src/atoms/molecule_test.cpp )
target_link_libraries( test_molecule libirc )
add_test( UNITTEST_molecule bin/test_molecule )

add_executable( test_io src/io/io_test.cpp )
target_link_libraries( test_io libirc )
add_test( UNITTEST_io bin/test_io )

add_executable( test_connectivity src/connectivity/connectivity_test.cpp )
target_link_libraries( test_connectivity libirc )
add_test( UNITTEST_connectivity bin/test_connectivity )

add_executable( test_wilson src/connectivity/wilson_test.cpp )
target_link_libraries( test_wilson libirc )
add_test( UNITTEST_wilson bin/test_wilson )

add_executable( test_transformation src/transformation/transformation_test.cpp )
target_link_libraries( test_transformation libirc )
add_test( UNITTEST_transformation bin/test_transformation )

add_executable( test_irc src/irc/irc_test.cpp )
target_link_libraries( test_irc libirc )
add_test( UNITTEST_irc bin/test_irc )