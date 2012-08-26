#-------------------------------------------------------------------------------
# serment/src/erme_geometry/CMakeLists.txt
#-------------------------------------------------------------------------------

SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(SRC
    Node.cc
    CartesianNode.cc
    CartesianNodeDetran.cc
    NodeFactoryDetran.cc
    NodePartitioner.cc
)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

ADD_LIBRARY(erme_geometry
            SHARED
            ${SRC} 
)
TARGET_LINK_LIBRARIES( erme_geometry
                       serment_comm
                       ${Detran_LIBRARIES}                    
)

INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/lib erme_geometry)

#-------------------------------------#
# TESTING
#-------------------------------------#

if( SERMENT_ENABLE_TEST )
    include_directories(${SRC_DIR}/test)
    add_subdirectory(${SRC_DIR}/test)
endif()

