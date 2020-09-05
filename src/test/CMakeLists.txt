# + ------------- +
# | CONFIGURATION |
# + ------------- +

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/config.h
)

# + ------- +
# | TESTING |
# + ------- +
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/)

set(TEST_EXECUTABLE_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/mathtools_test.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/linalg_test.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/periodic_table_test.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/atom_test.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molecule_test.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/io_test.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/connectivity_test.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/wilson_test.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/transformation_test.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/irc_test.cpp
)

foreach(TEST_SOURCE ${TEST_EXECUTABLE_SOURCES})
  get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
  add_executable(${TEST_NAME} ${TEST_SOURCE})
  target_link_libraries(${TEST_NAME} irc catch_defined)
  target_include_directories(${TEST_NAME} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()
