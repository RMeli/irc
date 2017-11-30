cmake_minimum_required(VERSION 3.8)
project(irc)

set(CMAKE_CXX_STANDARD 11)

add_library( libirc SHARED
                    src/atoms/atom.h
                    src/atoms/atom.cpp
                    src/atoms/molecule.h
                    src/atoms/molecule.cpp
                    src/atoms/periodic_table.h
                    src/atoms/periodic_table.cpp
                    src/connectivity/connectivity.h
                    src/connectivity/connectivity.cpp
        src/tools/comparison.h
        src/tools/comparison.cpp)
set_target_properties(libirc PROPERTIES OUTPUT_NAME "irc")
set_target_properties(libirc PROPERTIES LINKER_LANGUAGE CXX)


add_executable(irc src/main.cpp)
target_link_libraries( irc libirc )

find_package(Boost
        1.65.0
        REQUIRED
        )
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(libirc ${Boost_LIBRARIES})

find_package(Armadillo
        8.200.0
        REQUIRED
        )
include_directories(${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(libirc ${ARMADILLO_LIBRARIES})

enable_testing()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable( test_comparison src/tools/comparison_test.cpp)
target_link_libraries( test_comparison libirc )
add_test( UNITTEST_comparison bin/test_comparison )

add_executable( test_periodic_table src/atoms/periodic_table_test.cpp )
target_link_libraries( test_periodic_table libirc )
add_test( UNITTEST_periodic_table bin/test_periodic_table )

add_executable( test_atom src/atoms/atom_test.cpp )
target_link_libraries( test_atom libirc )
add_test( UNITTEST_atom bin/test_atom )

add_executable( test_molecule src/atoms/molecule_test.cpp )
target_link_libraries( test_molecule libirc )
add_test( UNITTEST_molecule bin/test_molecule )

add_executable( test_connectivity src/connectivity/connectivity_test.cpp )
target_link_libraries( test_connectivity libirc )
add_test( UNITTEST_connectivity bin/test_connectivity )