#-------------------------------------------------------------------------------
# serment/src/erme_solvers/CMakeLists.txt
#-------------------------------------------------------------------------------

SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(SRC
    GlobalSolverBase.cc
    GlobalSolverPicard.cc
    GlobalSolverNewton.cc
    GlobalSolverSteffensen.cc
    Jacobian.cc
    FullJacobian.cc
    ManagerERME.cc    
    NonlinearResidual.cc
)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

ADD_LIBRARY(erme_solver
            SHARED
            ${SRC} 
)
TARGET_LINK_LIBRARIES(erme_solver
                      erme             
)

INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/lib erme_solver)

#-------------------------------------#
# TESTING
#-------------------------------------#

if(SERMENT_ENABLE_TEST)
  include_directories(${SRC_DIR}/test)
  add_subdirectory(${SRC_DIR}/test)
endif()

