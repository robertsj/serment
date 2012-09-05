#-------------------------------------------------------------------------------
# serment/src/erme/CMakeLists.txt
#-------------------------------------------------------------------------------

SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(SRC
    StateERME.cc
    Connect.cc
    ResponseMatrix.cc
    FissionOperator.cc
    AbsorptionOperator.cc
    LeakageOperator.cc
)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

ADD_LIBRARY(erme
            SHARED
            ${SRC} 
)
TARGET_LINK_LIBRARIES( erme
                       erme_response
                       erme_geometry
                       linear_algebra
                       serment_comm
                       ${Detran_LIBRARIES}                    
)

INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/lib erme)

#-------------------------------------#
# TESTING
#-------------------------------------#

if( SERMENT_ENABLE_TEST )
    include_directories(${SRC_DIR}/test)
    add_subdirectory(${SRC_DIR}/test)
endif()

