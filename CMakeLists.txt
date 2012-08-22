#-------------------------------------------------------------------------------
# serment/CMakeLists.txt
#
# Note, this file was based in part on S. Johnson's PyTRT's CMake-based build
# system, which, like the rest of the code, is quite nicely written.
#-------------------------------------------------------------------------------


##===========================================================================##
## PROJECT                                                                   ##
##===========================================================================##

project( SERMENT )
cmake_minimum_required( VERSION 2.6 FATAL_ERROR )
set(SERMENT_DIR ${CMAKE_CURRENT_SOURCE_DIR})


##===========================================================================##
## VERSION                                                                   ##
##===========================================================================##

set( VERSION_MAJOR 0 )
set( VERSION_MINOR 1 )
set( VERSION_PATCH_LEVEL 0 )
message( 
"++ SERMENT version ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH_LEVEL}" 
)


##===========================================================================##
## CMAKE MODULES                                                             ##
##===========================================================================##

# CMake files included with SERMENT
LIST( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/tools/CMake" )


##===========================================================================##
## USER OPTIONS                                                              ##
##===========================================================================##

# Set the default build type (debug)
if( NOT CMAKE_BUILD_TYPE)
    set( CMAKE_BUILD_TYPE Debug CACHE STRING
        "Choose the type of build, options are: Debug Release." FORCE)
endif()

# Options for DBC macros.  For now, all on or all off.  Note, there *is* a
#   performance hit for having it on!
set( SERMENT_ENABLE_DBC YES CACHE BOOL "Enable design by contract macros.")
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    if(SERMENT_ENABLE_DBC)
        message("WARNING: using DBC macros with Release configuration!")
    endif()
endif()

# Other build options.
set( SERMENT_ENABLE_DEBUG     YES CACHE BOOL  "Enable diagonostic print outs."                   )
set( SERMENT_ENABLE_TEST      NO CACHE BOOL  "Enable testing."                                   )
set( SERMENT_BUILD_ERME       YES CACHE BOOL  "Build eigenvalue response matrix equation solver.")
set( SERMENT_BUILD_DIFF2D     YES CACHE BOOL  "Build two dimensional diffusion solver."          )

# Options for using external libraries.
set(SERMENT_ENABLE_PETSC      YES CACHE BOOL  "Use PETSc solvers and datatypes."                 )
set(SERMENT_ENABLE_SLEPC      YES CACHE BOOL  "Use SLEPc solvers and datatypes."                 )
set( SERMENT_ENABLE_SILO      YES CACHE BOOL  "Use the Silo library for data i/o."               )
set( SERMENT_ENABLE_XML       YES CACHE BOOL  "Use the libxml2 library for input."               )
set( SERMENT_ENABLE_MPI       NO CACHE BOOL   "Enable MPI. (not used)"                           )
set( SERMENT_ENABLE_OPENMP    NO CACHE BOOL   "Enable OpenMP. (not used)"                        )

# PETSc is essential for most things...
if( NOT SERMENT_ENABLE_PETSC )
  if( SERMENT_BUILD_ERME OR SERMENT_BUILD_DIFF2D )
    message( FATAL_ERROR "You must enable PETSc for ERME and DIFF2D." )
  endif()
endif()


##===========================================================================##
## EXTERNAL LIBRARIES                                                        ##
##===========================================================================##

# Detran
find_package(Detran REQUIRED)
include_directories(${Detran_INCLUDE_DIR})
message("Detran library = ${Detran_LIBRARIES}")

# PETSc
if(SERMENT_ENABLE_PETSC)
    find_package(PETSc REQUIRED)
    include_directories(${PETSC_INCLUDES})
    add_definitions(${PETSC_DEFINITIONS})
endif()

# SLEPc (optional)
if(SERMENT_ENABLE_SLEPC)
    find_package(SLEPc REQUIRED)
    if(HAVE_SLEPC)
        include_directories(${SLEPC_INCLUDES})
    else()
        set( SERMENT_ENABLE_SLEPC OFF )
        message( STATUS "Disabling SLEPc since it was not found.  Try setting
                 the environment variable SLEPC_DIR." )
    endif()
endif()

if(SERMENT_ENABLE_SILO)
    find_package(Silo)
    if(SILO_FOUND)
        include_directories(${Silo_INCLUDE_DIR})
    else()
        set( SERMENT_ENABLE_SILO OFF )
        MESSAGE( STATUS "Disabling Silo since it was not found." )
    endif()
endif()

if(SERMENT_ENABLE_XML)
    find_package(LibXml2 REQUIRED)
    include_directories(${LIBXML2_INCLUDE_DIR} )
    add_definitions(${LIBXML2_DEFINITIONS})
endif()

if(SERMENT_ENABLE_MPI)
  find_package(MPI REQUIRED)
  include_directories(${MPI_INCLUDE_PATH})
endif()

if(SERMENT_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()


##===========================================================================##
## CONFIGURATION                                                             ##
##===========================================================================##

# Write config file, install it, and include that directory for all source
set( SERMENT_CONFIG_IN  "${CMAKE_CURRENT_SOURCE_DIR}/serment_config.h.in" )
set( SERMENT_CONFIG_OUT "${PROJECT_BINARY_DIR}/config/serment_config.h" )
configure_file( ${SERMENT_CONFIG_IN} ${SERMENT_CONFIG_OUT} )
include_directories( "${PROJECT_BINARY_DIR}/config" )


##===========================================================================##
## TESTING                                                                   ##
##===========================================================================##

option( USE_TESTING "Build testing" OFF )
if ( SERMENT_ENABLE_TEST )
    SET ( CTEST_INITIAL_CACHE "
      COVERAGE_COMMAND:FILEPATH=/usr/bin/gcov
    " )
    # need for coverage -- probably a cleaner way to set this up.
    set( CMAKE_CXX_FLAGS_DEBUG        
         "-g -O0  -fprofile-arcs -ftest-coverage" CACHE STRING 
         "Choose"
         FORCE
    )
    set( CMAKE_CXX_FLAGS        
         "-g -O0  -fprofile-arcs -ftest-coverage" CACHE STRING
         "Choose"
         FORCE
    )
    set( CMAKE_C_FLAGS_DEBUG          
         "-g -O0  -fprofile-arcs -ftest-coverage" CACHE STRING
         "Choose"
         FORCE
    )
    set( CMAKE_EXE_LINKER_FLAGS 
         "-g -O0  -fprofile-arcs -ftest-coverage" CACHE STRING
         "Choose"
         FORCE
    )
    set( CTEST_MEMORYCHECK_COMMAND 
       /home/robertsj/opt/valgrind/3.6.1/bin/valgrind CACHE STRING 
       "Set the memory checking command."
       FORCE
    )
    set( CTEST_MEMORYCHECK_COMMAND_OPTIONS 
       "--tool=memcheck --leak-check=full" CACHE STRING 
       "Set the memory checking command options."
       FORCE
    ) 
    enable_testing()
    include(CTest)
    include_directories(src/linear_algebra/test)
endif()


##===========================================================================##
## SOURCE                                                                    ##
##===========================================================================##

include_directories( src )
add_subdirectory( src )


##===========================================================================##
## INSTALL                                                                    ##
##===========================================================================##

if( NOT CMAKE_INSTALL_PREFIX)
  set( CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR} CACHE STRING
       "Installation prefix; default is build directory."
       FORCE)
endif()

##===========================================================================##
## SOME DIAGNOSTIC OUTPUT                                                    ##
##===========================================================================##

message("++++ System:      ${CMAKE_SYSTEM} ${CMAKE_SYSTEM_PROCESSOR}")
message("++++ CXX:         ${CMAKE_CXX_COMPILER}")
message("++++ Build type:  ${CMAKE_BUILD_TYPE}")

if(SERMENT_ENABLE_DBC)
  message("++++++ DBC enabled")
else()
  message("++++++ DBC DISABLED")
endif()

if(SERMENT_ENABLE_PETSC)
  message("++++++ PETSc enabled")
else()
  message("++++++ PETSc DISABLED")
endif()
if(SERMENT_ENABLE_SLEPC)
  message("++++++ SLEPc enabled")
else()
  message("++++++ SLEPc DISABLED")
endif()
if(SERMENT_ENABLE_SILO)
  message("++++++ Silo enabled")
else()
  message("++++++ Silo DISABLED")
endif()
if(SERMENT_ENABLE_MPI)
  message("++++++ MPI enabled")
else()
  message("++++++ MPI DISABLED")
endif()
if(SERMENT_ENABLE_OPENMP)
  message("++++++ OpenMP enabled")
else()
  message("++++++ OpenMP DISABLED")
endif()

if(NOT SERMENT_BUILD_ERME)
  message("++++++ erme DISABLED")
endif()
if(NOT SERMENT_BUILD_DIFF2D)
  message("++++++ diff2d DISABLED")
endif()
if(NOT SERMENT_BUILD_SLABTRAN)
  message("++++++ slabtran DISABLED")
endif()



