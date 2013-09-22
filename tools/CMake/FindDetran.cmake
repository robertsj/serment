# This module defines
#   Detran_INCLUDE_DIR
#   Detran_LIBRARIES, libraries to link against to use Detran.
#   Detran_FOUND, If false, do not try to use Detran.
#
# The user should specify the head Detran directory, Detran_DIR.

if (Detran_DIR)
    message( "-- Using Detran_DIR: ${Detran_DIR} " )
endif()

# Set the cmake library path.  This is because a clean
# build of cmake 2.8.4 seems to exclude this path.  It might
# be a Ubuntu-specific quirk.
set(CMAKE_LIBRARY_PATH "/usr/lib/x86_64-linux-gnu")

find_path(Detran_INCLUDE_DIR 
          NAMES detran_config.hh
          PATHS ${Detran_DIR}/include 
                ${Detran_INC}
)
find_library(Detran_utilities   NAMES utilities   PATHS ${Detran_DIR}/lib ${Detran_LIB})
find_library(Detran_callow      NAMES callow      PATHS ${Detran_DIR}/lib ${Detran_LIB})
find_library(Detran_angle       NAMES angle       PATHS ${Detran_DIR}/lib ${Detran_LIB})
find_library(Detran_orthog      NAMES orthog      PATHS ${Detran_DIR}/lib ${Detran_LIB})
find_library(Detran_geometry    NAMES geometry    PATHS ${Detran_DIR}/lib ${Detran_LIB})
find_library(Detran_material    NAMES material    PATHS ${Detran_DIR}/lib ${Detran_LIB})
find_library(Detran_boundary    NAMES boundary    PATHS ${Detran_DIR}/lib ${Detran_LIB})
find_library(Detran_transport   NAMES transport   PATHS ${Detran_DIR}/lib ${Detran_LIB})
find_library(Detran_solvers     NAMES solvers     PATHS ${Detran_DIR}/lib ${Detran_LIB})
find_library(Detran_ioutils     NAMES ioutils     PATHS ${Detran_DIR}/lib ${Detran_LIB})
find_library(Detran_postprocess NAMES postprocess PATHS ${Detran_DIR}/lib ${Detran_LIB})
find_library(Detran_kinetics    NAMES kinetics PATHS ${Detran_DIR}/lib ${Detran_LIB})
find_library(Detran_external_source  NAMES external_source PATHS ${Detran_DIR}/lib ${Detran_LIB})


# handle the QUIETLY and REQUIRED arguments and set Detran_FOUND to TRUE if 
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Detran
  "Detran could not be found.  Try setting Detran_DIR or Detran_LIB and Detran_INC."
  Detran_INCLUDE_DIR 
  Detran_utilities
  Detran_callow
  Detran_angle
  Detran_orthog
  Detran_geometry
  Detran_material 
  Detran_boundary 
  Detran_transport 
  Detran_solvers
  Detran_ioutils
  Detran_postprocess
  Detran_kinetics
  Detran_external_source
)

if(DETRAN_FOUND)
  mark_as_advanced(Detran_INCLUDE_DIR Detran_LIBRARY)
  set(Detran_LIBRARIES ${Detran_utilities} 
                       ${Detran_callow} 
                       ${Detran_angle} 
                       ${Detran_orthog} 
                       ${Detran_geometry} 
                       ${Detran_material} 
                       ${Detran_boundary} 
                       ${Detran_transport} 
                       ${Detran_solvers}
                       ${Detran_ioutils}
                       ${Detran_postprocess}
                       ${Detran_kinetics}
                       ${Detran_external_source}
  )
  message("Detran IS FOUND")
else(DETRAN_FOUND)
  message("Detran NOT FOUND")
endif(DETRAN_FOUND)

message("DETRAN LIBRARIES ARE ${Detran_LIBRARIES}")


