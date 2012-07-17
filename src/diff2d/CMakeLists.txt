#-------------------------------------------------------------------------------
# serment/src/diff2d/CMakeLists.txt
#
# $Rev:: 103                                           $:Rev of last commit
# $Author:: j.alyn.roberts@gma                         $:Author of last commit
# $Date:: 2011-06-14 22:22:46 -0400 (Tue, 14 Jun 2011) $:Date of last commit  
#
#-------------------------------------------------------------------------------

SET( DIFF2D_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# Source subdirectories.
SET( DIFF2D_SRC 
    Diff2dInput.cc 
    Diff2dOutput.cc
    Diff2dProblem.cc
    Diff2dSolver.cc
)

# Include this directory.
INCLUDE_DIRECTORIES( ${DIFF2D_DIR} )

# Make library.
ADD_LIBRARY( diff2d
             ${DIFF2D_SRC} 
)
TARGET_LINK_LIBRARIES( diff2d
                       linalg
                       utils
                       ${PETSC_LIBRARIES} 
                       ${Silo_LIBRARIES}
)

# Make executable if desired
IF(SERMENT_BUILD_DIFF2D)
    ADD_EXECUTABLE( diff2dx diff2d.cc )
    TARGET_LINK_LIBRARIES( diff2dx
                           diff2d
                           linalg
                           utils
                           ${PETSC_LIBRARIES} 
                           ${Silo_LIBRARIES}
    )
ENDIF()

