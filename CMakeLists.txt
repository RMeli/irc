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

option(COVERAGE "Use code coverage if compiler is Clang" OFF)

if(COVERAGE)
    set(COVERAGE_COMPILE_FLAGS "--coverage -g -O0 -fno-inline -fno-inline-small-functions -fno-default-inline")
    set(COVERAGE_LINK_FLAGS    "--coverage")

    message("-- Code coverage enabled: ${COVERAGE_COMPILE_FLAGS}")

    set( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${COVERAGE_COMPILE_FLAGS}" )
    set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${COVERAGE_COMPILE_FLAGS}" )
    set( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${COVERAGE_LINK_FLAGS}" )
endif()

# + ----------- +
# | DIRECTORIES |
# + ----------- +

set(IRC_MOLECULES_DIR ${CMAKE_SOURCE_DIR}/molecules/)

# + ----------- +
# | IRC LIBRARY |
# + ----------- +

add_library( irc INTERFACE )

# + -------- +
# | INCLUDES |
# + -------- +

target_include_directories( irc INTERFACE include )
target_include_directories( irc INTERFACE include/libirc )
target_include_directories( irc INTERFACE ${CMAKE_BINARY_DIR}/src)

target_include_directories( irc INTERFACE external )

# + ----------------------- +
# | FIND EXTERNAL LIBRARIES |
# + ----------------------- +

#if(NOT EXISTS ${CMAKE_SOURCE_DIR}/include/boost)
find_package(
        Boost
        1.58.0
        REQUIRED
    )
include_directories(irc ${Boost_INCLUDE_DIRS}) # DON'T USE
target_link_libraries(irc INTERFACE ${Boost_LIBRARIES})
#endif()

option(WITH_ARMA "Use Armadillo for the linear algebra library" OFF)
option(WITH_EIGEN "Use Eigen for the linear algebra library" OFF)

if(WITH_ARMA)
    find_package(
        Armadillo
        6.500.0
        REQUIRED
    )

    if( ${ARMADILLO_FOUND})
        include_directories(irc ${ARMADILLO_INCLUDE_DIRS}) # DON'T USE
        target_link_libraries(irc INTERFACE ${ARMADILLO_LIBRARIES})
        add_definitions(-DHAVE_ARMA)

        message("-- Armadillo version: " ${ARMADILLO_VERSION_STRING})
    endif()
elseif(WITH_EIGEN)
    find_package(
        Eigen3
        3.2
        REQUIRED
    )

    if( ${EIGEN3_FOUND} )
        include_directories(${EIGEN3_INCLUDE_DIRS})
        target_link_libraries(irc INTERFACE ${EIGEN3_LIBRARIES})
        add_definitions(-DHAVE_EIGEN3)

        add_definitions(
                -DEIGEN_MATRIX_PLUGIN="eigen/plugins/Matrix_initializer_list.h"
        )

        message("-- Eigen3 version: " ${EIGEN3_VERSION_STRING})
    endif()
else()
    message( FATAL_ERROR "No linera algebra library provided. Specifiy -DWITH_ARMA or -DWITH_EIGEN" )
endif()

# + ------- +
# | TESTING |
# + ------- +

enable_testing()
add_subdirectory(src/test)
