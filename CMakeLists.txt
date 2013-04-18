cmake_minimum_required(VERSION 2.8)

project(ODFfibConvert)

find_package(ITK REQUIRED)
IF(ITK_FOUND)
  include(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  message(FATAL_ERROR "ITK not found. Please set ITK_DIR")
ENDIF(ITK_FOUND)

find_package(GenerateCLP REQUIRED)
IF(GenerateCLP_FOUND)
  include(${GenerateCLP_USE_FILE})
ELSE(GenerateCLP_FOUND)
  message(FATAL_ERROR "GenerateCLP not found. Please set GenerateCLP_DIR")
ENDIF(GenerateCLP_FOUND)

configure_file( ODFfibConvert.cpp.in ${CMAKE_CURRENT_BINARY_DIR}/ODFfibConvert.cpp )  # for the location of datFile & FacesFile
include_directories( ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ) # So the main cpp finds the headers
set( sources ${CMAKE_CURRENT_BINARY_DIR}/ODFfibConvert.cpp zconf.h zlib.h ODFReconstructor.h ODFReconstructor.cxx )

GENERATECLP(sources ODFfibConvert.xml)

add_executable(ODFfibConvert ${sources})

target_link_libraries(ODFfibConvert ${ITK_LIBRARIES} z) # add z library (zlib.h) to open gz files
