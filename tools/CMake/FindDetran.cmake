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
          NAMES DBC.hh
          PATHS ${Detran_DIR}/include 
                ${Detran_INC}
)

find_library(Detran_LIBRARY
             NAMES utilities 
             PATHS ${Detran_DIR}/lib
                   ${Detran_LIB}
)

# handle the QUIETLY and REQUIRED arguments and set Detran_FOUND to TRUE if 
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Detran 
  "Detran could not be found.  Try setting Detran_DIR or Detran_LIB and Detran_INC."
  Detran_LIBRARY Detran_INCLUDE_DIR)

if(DETRAN_FOUND)
  mark_as_advanced(Detran_INCLUDE_DIR Detran_LIBRARY)
  set(Detran_LIBRARIES ${Detran_LIBRARY})
  message("Detran IS FOUND")
else(DETRAN_FOUND)
message("Detran NOT FOUND")
endif(DETRAN_FOUND)


