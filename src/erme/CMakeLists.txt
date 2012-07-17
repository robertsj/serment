#-------------------------------------------------------------------------------
# serment/src/erme/CMakeLists.txt
#
# $Rev:: 103                                           $:Rev of last commit
# $Author:: j.alyn.roberts@gma                         $:Author of last commit
# $Date:: 2011-06-14 22:22:46 -0400 (Tue, 14 Jun 2011) $:Date of last commit  
#
#-------------------------------------------------------------------------------

SET( ERME_DIR ${CMAKE_CURRENT_SOURCE_DIR})

SET( ERME_SRC
# driver
    driver/GlobalInput.cc
# geometry
    geometry/Connect2dCart.cc
# global
#
# local
    local/LocalProblemDiff2d.cc
# response
    response/AbsorptionResponse.cc
    response/FissionResponse.cc
    response/LeakageResponse.cc
    response/ResponseFunction.cc
    response/ResponseFunctionDiffusion.cc
    response/ResponseFunctionServer.cc
    response/ResponseMatrix.cc
    response/ResponseMatrixFull.cc
# solvers
    solvers/GlobalSolver.cc
    solvers/InnerIterPower.cc
    solvers/InnerIterSLEPc.cc
    solvers/InvItShell.cc
    solvers/JacobianShell.cc
    solvers/Newton.cc
    solvers/PCAppxJacobian.cc
    solvers/PowerIter.cc
)

# Add include directories
INCLUDE_DIRECTORIES( ${ERME_DIR}/driver )
INCLUDE_DIRECTORIES( ${ERME_DIR}/geometry )
INCLUDE_DIRECTORIES( ${ERME_DIR}/global )
INCLUDE_DIRECTORIES( ${ERME_DIR}/local )
INCLUDE_DIRECTORIES( ${ERME_DIR}/response )
INCLUDE_DIRECTORIES( ${ERME_DIR}/solvers )

# Make the library
ADD_LIBRARY( erme
             ${ERME_SRC} 
)
TARGET_LINK_LIBRARIES( erme
                       diff2d
                       utils
                       ${PETSC_LIBRARIES} 
                       ${SLEPC_LIBRARIES}
                       ${SILO_LIBRARIES} 
)

# Build the serment executable
if( SERMENT_BUILD_ERME )
    ADD_EXECUTABLE( sermentx driver/Serment.cc )
    TARGET_LINK_LIBRARIES( sermentx 
                           erme
                           diff2d
                           utils
                           ${PETSC_LIBRARIES} 
                           ${SLEPC_LIBRARIES}
                           ${SILO_LIBRARIES} 
    )
endif()

INSTALL_TARGETS(/bin sermentx)
