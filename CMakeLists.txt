cmake_minimum_required(VERSION 2.8)

project(ODFfibConvert)

# ITK
find_package(ITK REQUIRED)
IF(ITK_FOUND)
  include(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  message(FATAL_ERROR "ITK not found. Please set ITK_DIR")
ENDIF(ITK_FOUND)

# GenerateCLP
find_package(GenerateCLP REQUIRED)
IF(GenerateCLP_FOUND)
  include(${GenerateCLP_USE_FILE})
ELSE(GenerateCLP_FOUND)
  message(FATAL_ERROR "GenerateCLP not found. Please set GenerateCLP_DIR")
ENDIF(GenerateCLP_FOUND)

# Qt4
find_package(Qt4 REQUIRED) # Qt needed to compile headers of DSIstudio
if(QT_USE_FILE)
  include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${QT_INCLUDE_DIR})
  include(${QT_USE_FILE})
  add_definitions(-DQT_GUI_LIBS -DQT_CORE_LIB -DQT3_SUPPORT)
else(QT_USE_FILE)
  message(FATAL_ERROR, "QT not found. Please set QT_DIR.")
endif(QT_USE_FILE)

# CLAPACK for ODF reconstruction (see Spharm-PDM)
find_package(CLAPACK NO_MODULE REQUIRED)
include_directories(${CLAPACK_DIR}/../CLAPACK/INCLUDE) # ${CLAPACK_DIR} is the binary dir -> include files are in the source code
set(CLAPACK_LIBRARIES lapack blas f2c )

# DSI studio & Boost headers
set( DSIstudio_PATH "" CACHE PATH "Path to the DSIstudio source tree")
set( Boost_INCLUDE "" CACHE PATH "Path to the Boost Include directory") # /boost/ptr_container/ptr_vector.hpp
set( LibsFlags "-I${DSIstudio_PATH} -I${Boost_INCLUDE}" )
set( CMAKE_CXX_FLAGS    "${CMAKE_CXX_FLAGS} ${LibsFlags}" )
set( CMAKE_C_FLAGS      "${CMAKE_C_FLAGS} ${LibsFlags}" )
set( CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} ${LibsFlags}" )

# Sources
configure_file( define.h.in ${CMAKE_CURRENT_BINARY_DIR}/define.h )  # for the location of datFile & FacesFile
include_directories( ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ) # So the main cxx finds the headers
set( sources ODFfibConvert.cxx ODFSpharmComputer.cxx ODFReconstructor.cxx )
GENERATECLP(sources ODFfibConvert.xml)

# Executable
add_executable(ODFfibConvert ${sources})
target_link_libraries(ODFfibConvert ${QT_LIBRARIES} ${ITK_LIBRARIES} z ${CLAPACK_LIBRARIES}) # add z library (zlib.h) to open gz files
