cmake_minimum_required(VERSION 2.8)

project(simu)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

set(CMAKE_BUILD_TYPE "Release")
message("-- Build type: ${CMAKE_BUILD_TYPE}")

set(Boost_ADDITIONAL_VERSIONS "1.57" "1.57.0" "1.56" "1.56.0" "1.55" "1.55.0" "1.54" "1.54.0" "1.53" "1.53.0" "1.52" "1.52.0" "1.51" "1.51.0" "1.50" "1.50.0" "1.49" "1.49.0" "1.48" "1.48.0" "1.47" "1.47.0" "1.46" "1.46.0" "1.45" "1.45.0" "1.44" "1.44.0" "1.42" "1.42.0" "1.41.0" "1.41" "1.40.0" "1.40")

# find boost
find_package(Boost REQUIRED)
if (NOT Boost_FOUND)
  message(SEND_ERROR "Failed to find boost libraries.")
endif (NOT Boost_FOUND)

set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME_LIBS   ON)
set(Boost_USE_STATIC ON)
set(Boost_USE_STATIC_LIBS   ON)
set(SIMU_BOOST_COMPONENTS program_options filesystem system)
find_package(Boost COMPONENTS REQUIRED ${SIMU_BOOST_COMPONENTS})
if (NOT Boost_FOUND)
  message(SEND_ERROR "Failed to find required boost libraries.")
  return()
else (NOT Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif (NOT Boost_FOUND)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if (COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif (COMPILER_SUPPORTS_CXX0X)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else (COMPILER_SUPPORTS_CXX11)
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif (COMPILER_SUPPORTS_CXX11)

SET(DISABLE_TESTS ON CACHE BOOL "Disables building tests.")
SET(DISABLE_SHARED_LIBRARY OFF CACHE BOOL "Disables building shared libraries.")

add_subdirectory( libplinkio )
add_subdirectory( src )
