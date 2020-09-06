cmake_minimum_required(VERSION 3.3)
project(libirc VERSION 0.1.0 LANGUAGES CXX)

# + ----------------- +
# | BUILDING SETTINGS |
# + ----------------- +

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# + ----------------- +
# | LANGUAGE SETTINGS |
# + ----------------- +

# C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# + ----------------- +
# | COMPILATION FLAGS |
# + ----------------- +

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wall WALL_SUPPORTED)
if(${WALL_SUPPORTED})
  add_compile_options(-Wall)
endif()

option(COVERAGE "Use code coverage if compiler is Clang" OFF)
if(COVERAGE)
  set(COVERAGE_COMPILE_FLAGS "--coverage -g -O0 -fno-inline -fno-inline-small-functions -fno-default-inline")
  set(COVERAGE_LINK_FLAGS "--coverage")

  message("-- Code coverage enabled: ${COVERAGE_COMPILE_FLAGS}")

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COVERAGE_COMPILE_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COVERAGE_COMPILE_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${COVERAGE_LINK_FLAGS}")
endif()

# + ----------- +
# | DIRECTORIES |
# + ----------- +

set(IRC_MOLECULES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/molecules/)

# + ----------- +
# | IRC LIBRARY |
# + ----------- +

add_library(irc INTERFACE)

# + -------- +
# | INCLUDES |
# + -------- +

target_include_directories(irc INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

# + ----------------------- +
# | FIND EXTERNAL LIBRARIES |
# + ----------------------- +

find_package(Boost 1.58.0 REQUIRED)
if(TARGET Boost::headers)
  target_link_libraries(irc INTERFACE Boost::headers)
elseif(TARGET Boost::boost)
  target_link_libraries(irc INTERFACE Boost::boost)
else()
  message(WARNING "Boost inclusion path is hardcoded. Any installed package will be non-relocatable! Update to a newer boost and/or CMake version to fix this.")
  target_include_directories(irc INTERFACE ${Boost_INCLUDE_DIRS})
endif()

option(WITH_ARMA "Use Armadillo for the linear algebra library" OFF)
option(WITH_EIGEN "Use Eigen for the linear algebra library" OFF)

if(WITH_ARMA)
  find_package(Armadillo 6.500.0 REQUIRED)

  if(TARGET armadillo)
    target_link_libraries(irc INTERFACE armadillo)
  else()
    message(WARNING "Armadillo inclusion path is hardcoded. Any installation will be non-relocatable! Update to a newer CMake version to fix this.")
    target_include_directories(irc INTERFACE ${ARMADILLO_INCLUDE_DIRS})
    target_link_libraries(irc INTERFACE ${ARMADILLO_LIBRARIES})
  endif()

  target_compile_definitions(irc INTERFACE -DHAVE_ARMA)
  message("-- Armadillo version: " ${ARMADILLO_VERSION_STRING})
elseif(WITH_EIGEN)
  find_package(Eigen3 3.0 REQUIRED)
  target_link_libraries(irc INTERFACE Eigen3::Eigen)
  target_compile_definitions(irc INTERFACE
    -DHAVE_EIGEN3
    -DEIGEN_MATRIX_PLUGIN="libirc/plugins/eigen/Matrix_initializer_list.h"
  )

  message("-- Eigen3 version: " ${EIGEN3_VERSION_STRING})
else()
  message(FATAL_ERROR "No linear algebra library selected. Specify -DWITH_ARMA or -DWITH_EIGEN")
endif()

# + ------- +
# | TESTING |
# + ------- +

enable_testing()
add_subdirectory(external)
add_subdirectory(src/test)

# + ------------ +
# | INSTALLATION |
# + ------------ +

# Copy headers folder
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/libirc
  DESTINATION include
)
# Create an export set
install(TARGETS irc EXPORT ircTargets)
include(CMakePackageConfigHelpers)
# Version file
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/irc-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)
# Config file
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/irc-config.cmake
  INSTALL_DESTINATION lib/cmake/irc
)
# Targets files
export(
  EXPORT ircTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/irc-targets.cmake
)
install(
  EXPORT ircTargets
  FILE irc-targets.cmake
  DESTINATION lib/cmake/irc
)
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/irc-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/irc-config-version.cmake
  DESTINATION lib/cmake/irc
)
