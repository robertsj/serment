#-------------------------------------------------------------------------------
# serment/src/linalg/CMakeLists.txt
#
# $Rev:: 103                                           $:Rev of last commit
# $Author:: j.alyn.roberts@gma                         $:Author of last commit
# $Date:: 2011-06-14 22:22:46 -0400 (Tue, 14 Jun 2011) $:Date of last commit  
#
#-------------------------------------------------------------------------------

SET( LINALG_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# Source subdirectories.
SET( LINALG_SRC 
    SermentMatrixBCRS.cc  
    SermentMatrixCRS.cc  
    SermentMatrixFull.cc  
    SermentMatrixShell.cc  
    SermentVector.cc
)

# Include this directory.
INCLUDE_DIRECTORIES( ${LINALG_DIR} )

# Make library.
ADD_LIBRARY( linalg 
             ${LINALG_SRC} 
)
TARGET_LINK_LIBRARIES( linalg
                       ${PETSC_LIBRARIES} 
                       ${SILO_LIBRARIES}
)

## Install library and headers.
#INSTALL( TARGETS
#         ${HEADERS}
#         DESTINATION
#         ${CMAKE_INSTALL_PREFIX}/include/diff2d
#)
#INSTALL( LIBRARY
#         ${SUBPROJECT_NAME}lib
#         DESTINATION
#         ${CMAKE_INSTALL_PREFIX}/lib/diff2d
#)


