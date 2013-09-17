#-------------------------------------------------------------------------------
# Serment CMake Configuration
#-------------------------------------------------------------------------------

#------------------------------------------------------------------------------#
# PROJECT
#------------------------------------------------------------------------------#

project(Serment)
cmake_minimum_required( VERSION 2.6 FATAL_ERROR )
set(SERMENT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#set (CMAKE_CXX_COMPILER /home/robertsj/opt/tau/x86_64/bin/tau_cxx.sh)
#set (CMAKE_C_COMPILER tau_cc.sh)

#------------------------------------------------------------------------------#
# VERSION
#------------------------------------------------------------------------------#

set( VERSION_MAJOR 1 )
set( VERSION_MINOR 0 )
set( VERSION_PATCH_LEVEL 0 )
message( 
"++ SERMENT version ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH_LEVEL}" 
)

#------------------------------------------------------------------------------#
# CMAKE MODULES
#------------------------------------------------------------------------------#

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/tools/CMake")

#------------------------------------------------------------------------------#
# USER OPTIONS
#------------------------------------------------------------------------------#

# Set the default build type (debug)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING
        "Choose the type of build, options are: Debug Release." FORCE)
endif()

set(SERMENT_ENABLE_TEST    NO CACHE BOOL "Enable testing.")
set(SERMENT_ENABLE_HDF5    NO CACHE BOOL "Enable HDF5.")
set(SERMENT_ENABLE_MPI     NO CACHE BOOL "Enable MPI.")
set(SERMENT_ENABLE_PYTHON  NO CACHE BOOL "Enable Python bindings via SWIG.")

#------------------------------------------------------------------------------#
# EXTERNAL LIBRARIES
#------------------------------------------------------------------------------#

# Detran (required)
find_package(Detran REQUIRED)
include_directories(${Detran_INCLUDE_DIR})
list(APPEND CMAKE_MODULE_PATH "${Detran_INCLUDE_DIR}/../tools/CMake")

# PETSc (required)
find_package(PETSc REQUIRED)
include_directories(${PETSC_INCLUDES})
add_definitions(${PETSC_DEFINITIONS})
message("PETSc library = ${PETSC_LIBRARIES}")
find_package(MPI REQUIRED)

# SLEPc (required)
find_package(SLEPc REQUIRED)
include_directories(${SLEPC_INCLUDES})
message("SLEPc library = ${SLEPC_LIBRARIES}")

# Boost (required)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREAD OFF)
set(Boost_ADDITIONAL_VERSIONS "1.50" "1.50.0" "1.51" "1.51.0")
find_package(Boost 1.50 COMPONENTS serialization program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# HDF5 (optional)
if(SERMENT_ENABLE_HDF5)
  find_package(HDF5)
  if(HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIR})
  else()
    set(SERMENT_ENABLE_HDF5 OFF)
    message(STATUS "Disabling HDF5 since it was not found.")
    # Need to reset them, since FindHDF5 fills them with "not founds", 
    # which leads to errors downstream.
    set(HDF5_LIBRARIES "")
    set(HDF_INCLUDE_DIR "")
  endif()
endif()

# MPI (optional)
if(SERMENT_ENABLE_MPI)
  find_package(MPI REQUIRED)
  if(MPI_FOUND)
    include_directories(${MPI_INCLUDE_PATH})
  else()
    set(SERMENT_ENABLE_MPI OFF)
    message(STATUS "Disabling MPI since it was not found.")
    # Need to reset them, since FindHDF5 fills them with "not founds", 
    # which leads to errors downstream.
    set(MPI_LIBRARIES "")
    set(MPI_INCLUDE_PATH "")
  endif()
endif()

if(SERMENT_ENABLE_PYTHON)
  find_package(SWIG REQUIRED)
  include(${SWIG_USE_FILE})
  find_package(PythonLibs REQUIRED)
  include_directories(${PYTHON_INCLUDE_PATH})
  include_directories(${PROJECT_BINARY_DIR}/config)
  #set(CMAKE_SWIG_FLAGS "-verbose")
  set(CMAKE_SWIG_FLAGS "-v")
endif()

#set(GPERFTOOLS_LIBRARIES "/home/robertsj/opt/gperftools/lib/libprofiler.so")
#set(GPERFTOOLS_INCLUDE   "/home/robertsj/opt/gperftools/include")
#include_directories(${GPERFTOOLS_INCLUDE})

#------------------------------------------------------------------------------#
# CONFIGURATION
#------------------------------------------------------------------------------#

include(GetTheDate)
get_the_date(SERMENT_COMPILED_D SERMENT_COMPILED_M SERMENT_COMPILED_Y)

set(SERMENT_WC_REVISION_HASH "unknown")
find_package(Git)
if(GIT_FOUND)
  GIT_WC_INFO(${CMAKE_CURRENT_SOURCE_DIR} "SERMENT")
endif()

# Write config file, install it, and include that directory for all source
set(SERMENT_CONFIG_IN  "${CMAKE_CURRENT_SOURCE_DIR}/serment_config.h.in")
set(SERMENT_CONFIG_OUT "${PROJECT_BINARY_DIR}/config/serment_config.h")
configure_file(${SERMENT_CONFIG_IN} ${SERMENT_CONFIG_OUT})
include_directories("${PROJECT_BINARY_DIR}/config")


#------------------------------------------------------------------------------#
# TESTING
#------------------------------------------------------------------------------#

option( USE_TESTING "Build testing" OFF )
if ( SERMENT_ENABLE_TEST )
    enable_testing()
    include(CTest)
    include_directories(src/linear_algebra/test)
endif()

#------------------------------------------------------------------------------#
# SOURCE
#------------------------------------------------------------------------------#

include_directories( src )
add_subdirectory( src )

#------------------------------------------------------------------------------#
# INSTALL
#------------------------------------------------------------------------------#

if( NOT CMAKE_INSTALL_PREFIX)
  set( CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR} CACHE STRING
       "Installation prefix; default is build directory."
       FORCE)
endif()

#------------------------------------------------------------------------------#
# OUTPUT
#------------------------------------------------------------------------------#

message(" ")
message("---------------------------------------------------------------")
message("              SUMMARY OF CMAKE CONFIGURATION")
message("---------------------------------------------------------------")
message("++++ System:      ${CMAKE_SYSTEM} ${CMAKE_SYSTEM_PROCESSOR}")
message("++++ CXX:         ${CMAKE_CXX_COMPILER}")
message("++++ Build type:  ${CMAKE_BUILD_TYPE}")
message("++++ Install:     ${CMAKE_INSTALL_PREFIX}")
message("++++ Git SHA1:    ${SERMENT_WC_REVISION}")
message(" ")
message("---------------------------------------------------------------")
message("                      BUILD OPTIONS")
message("---------------------------------------------------------------")
if(SERMENT_ENABLE_DBC)
  message("++++++ DBC:    enabled (because DBC in Detran is enabled)")
else()
  message("++++++ DBC:    disabled (because DBC in Detran is disabled)")
endif()
if(SERMENT_ENABLE_MPI)
  message("++++++ MPI:    enabled")
else()
  message("++++++ MPI:    disabled")
endif()
message("---------------------------------------------------------------")

