# + ------------- +
# | CATCH LIBRARY |
# + ------------- +
add_library(catch INTERFACE)
target_include_directories(catch INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/catch)

add_library(catch_defined catch_main.cpp)
target_link_libraries(catch_defined catch)
