cmake_minimum_required(VERSION 3.8)
project(irc)

set(CMAKE_CXX_STANDARD 11)

add_library( libirc SHARED
                    )
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