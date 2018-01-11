# + ------- +
# | TESTING |
# + ------- +

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# math
add_executable( test_mathtools mathtools_test.cpp test.cpp)
target_link_libraries( test_mathtools irc )
add_test( UNITTEST_mathtools ${CMAKE_BINARY_DIR}/bin/test_mathtools )

# periodic_table
add_executable( test_periodic_table periodic_table_test.cpp test.cpp)
target_link_libraries( test_periodic_table irc )
add_test( UNITTEST_periodic_table ${CMAKE_BINARY_DIR}/bin/test_periodic_table)

# atom
add_executable( test_atom atom_test.cpp test.cpp)
target_link_libraries( test_atom irc )
add_test( UNITTEST_atom ${CMAKE_BINARY_DIR}/bin/test_atom )

# molecule
add_executable( test_molecule molecule_test.cpp test.cpp)
target_link_libraries( test_molecule irc )
add_test( UNITTEST_molecule ${CMAKE_BINARY_DIR}/bin/test_molecule )

# io
add_executable( test_io io_test.cpp test.cpp)
target_link_libraries( test_io irc )
add_test( UNITTEST_io ${CMAKE_BINARY_DIR}/bin/test_io )

# connectivity
add_executable( test_connectivity connectivity_test.cpp test.cpp)
target_link_libraries( test_connectivity irc )
add_test( UNITTEST_connectivity ${CMAKE_BINARY_DIR}/bin/test_connectivity )

# wilson
add_executable( test_wilson wilson_test.cpp test.cpp)
target_link_libraries( test_wilson irc )
add_test( UNITTEST_wilson ${CMAKE_BINARY_DIR}/bin/test_wilson )

# transformation
add_executable( test_transformation transformation_test.cpp test.cpp)
target_link_libraries( test_transformation irc )
add_test( UNITTEST_transformation ${CMAKE_BINARY_DIR}/bin/test_transformation )

# irc
add_executable( test_irc irc_test.cpp test.cpp)
target_link_libraries( test_irc irc )
add_test( UNITTEST_irc ${CMAKE_BINARY_DIR}/bin/test_irc )