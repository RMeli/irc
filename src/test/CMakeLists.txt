# + ------------- +
# | CATCH LIBRARY |
# + ------------- +
add_library(catch INTERFACE)
target_include_directories(catch INTERFACE ${CMAKE_SOURCE_DIR}/external/catch/)

add_library(catch_defined test.cpp)
target_link_libraries(catch_defined catch)

# + ------- +
# | TESTING |
# + ------- +
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/)

# math
add_executable( test_mathtools mathtools_test.cpp)
target_link_libraries( test_mathtools irc catch_defined)
add_test( UNITTEST_mathtools ${CMAKE_BINARY_DIR}/bin/test_mathtools )

# periodic_table
add_executable( test_periodic_table periodic_table_test.cpp)
target_link_libraries( test_periodic_table irc catch_defined)
add_test( UNITTEST_periodic_table ${CMAKE_BINARY_DIR}/bin/test_periodic_table)

# atom
add_executable( test_atom atom_test.cpp)
target_link_libraries( test_atom irc catch_defined)
add_test( UNITTEST_atom ${CMAKE_BINARY_DIR}/bin/test_atom )

# molecule
add_executable( test_molecule molecule_test.cpp)
target_link_libraries( test_molecule irc catch_defined)
add_test( UNITTEST_molecule ${CMAKE_BINARY_DIR}/bin/test_molecule )

# io
add_executable( test_io io_test.cpp)
target_link_libraries( test_io irc catch_defined)
add_test( UNITTEST_io ${CMAKE_BINARY_DIR}/bin/test_io )

# connectivity
add_executable( test_connectivity connectivity_test.cpp)
target_link_libraries( test_connectivity irc catch_defined)
add_test( UNITTEST_connectivity ${CMAKE_BINARY_DIR}/bin/test_connectivity )

# wilson
add_executable( test_wilson wilson_test.cpp)
target_link_libraries( test_wilson irc catch_defined)
add_test( UNITTEST_wilson ${CMAKE_BINARY_DIR}/bin/test_wilson )

# transformation
add_executable( test_transformation transformation_test.cpp)
target_link_libraries( test_transformation irc catch_defined)
add_test( UNITTEST_transformation ${CMAKE_BINARY_DIR}/bin/test_transformation )

# irc
add_executable( test_irc irc_test.cpp)
target_link_libraries( test_irc irc catch_defined)
add_test( UNITTEST_irc ${CMAKE_BINARY_DIR}/bin/test_irc )
