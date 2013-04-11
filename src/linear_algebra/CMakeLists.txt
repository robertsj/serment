#-------------------------------------------------------------------------------
# serment/src/linear_algebra/CMakeLists.txt
#-------------------------------------------------------------------------------

SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(SRC
    MatrixBase.cc
    Matrix.cc
    MatrixShell.cc
    Vector.cc
    EigenSolver.cc
)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

ADD_LIBRARY(linear_algebra
            SHARED
            ${SRC} 
)
MESSAGE(" petsc =  ${PETSC_LIBRARIES}")
TARGET_LINK_LIBRARIES(linear_algebra
                      ${Detran_LIBRARIES}
                      ${PETSC_LIBRARIES}
                      ${SLEPC_LIBRARIES}
                      ${Boost_LIBRARIES}
)

INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/lib linear_algebra)

#-------------------------------------#
# TESTING
#-------------------------------------#

if( SERMENT_ENABLE_TEST )
    include_directories(${SRC_DIR}/test)
    add_subdirectory(${SRC_DIR}/test)
endif()

