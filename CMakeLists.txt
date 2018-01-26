cmake_minimum_required(VERSION 3.3)
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

# + ----------------- +
# | COMPILATION FLAGS |
# + ----------------- +

add_compile_options(-Wall) # DON'T USE

# + ----------- +
# | DIRECTORIES |
# + ----------- +

set(IRC_MOLECULES_DIR ${CMAKE_SOURCE_DIR}/molecules/)

# + ----------- +
# | IRC LIBRARY |
# + ----------- +

configure_file( ${CMAKE_SOURCE_DIR}/src/config.h.in
                ${CMAKE_BINARY_DIR}/src/config.h)

set(IRC_HEADERS
        include/libirc.h
        include/libirc/atom.h
        include/libirc/connectivity.h
        include/libirc/constants.h
        include/libirc/conversion.h
        include/libirc/io.h
        include/libirc/irc.h
        include/libirc/linalg.h
        include/libirc/mathtools.h
        include/libirc/molecule.h
        include/libirc/periodic_table.h
        include/libirc/transformation.h
        include/libirc/wilson.h)

set(IRC_SOURCES
        src/atoms/atom.cpp
        src/atoms/periodic_table.cpp
        src/tools/mathtools.cpp)

add_library( irc SHARED ${IRC_HEADERS} ${IRC_SOURCES})

set_target_properties(irc PROPERTIES LINKER_LANGUAGE CXX)

# + -------- +
# | INCLUDES |
# + -------- +

target_include_directories( irc PUBLIC include )
target_include_directories( irc PUBLIC include/libirc )
target_include_directories( irc PUBLIC ${CMAKE_BINARY_DIR}/src)

# + ----------------------- +
# | FIND EXTERNAL LIBRARIES |
# + ----------------------- +

if(NOT EXISTS ${CMAKE_SOURCE_DIR}/include/boost)
    find_package(
            Boost
            1.60.0
            REQUIRED
        )
    include_directories(irc ${Boost_INCLUDE_DIRS}) # DON'T USE
    target_link_libraries(irc ${Boost_LIBRARIES})
endif()

if(WITH_ARMA)
    find_package(
        Armadillo
        8.000.0
        REQUIRED
    )

    if( ${ARMADILLO_FOUND})
        include_directories(irc ${ARMADILLO_INCLUDE_DIRS}) # DON'T USE
        target_link_libraries(irc ${ARMADILLO_LIBRARIES})
        add_definitions(-DHAVE_ARMA)

        message("-- Armadillo version: " ${ARMADILLO_VERSION_STRING})
    endif()
elseif(WITH_EIGEN)
    find_package(
        Eigen3
        3.2
        REQUIRED
    )

    if( ${EIGEN3_FOUND})
        include_directories(${EIGEN3_INCLUDE_DIRS})
        target_link_libraries(irc ${EIGEN3_LIBRARIES})
        add_definitions(-DHAVE_EIGEN3)

        add_definitions(
                -DEIGEN_MATRIX_PLUGIN="eigen/plugins/Matrix_initializer_list.h"
        )
        add_definitions(
                -DEIGEN_ARRAY_PLUGIN="eigen/plugins/Array_initializer_list.h"
        )

        message("-- Eigen3 version: " ${EIGEN3_VERSION_STRING})
    endif()
else()
    message( FATAL_ERROR "No linera algebra library provided." )
endif()

# + ------- +
# | INSTALL |
# + ------- +

install (TARGETS irc
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib)

install(DIRECTORY include/
        DESTINATION include/${PROJECT_NAME}
        FILES_MATCHING PATTERN "*.h*")

# + ------- +
# | TESTING |
# + ------- +

enable_testing()
add_subdirectory(src/test)