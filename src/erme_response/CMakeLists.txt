#-------------------------------------------------------------------------------
# serment/src/erme_response/CMakeLists.txt
#-------------------------------------------------------------------------------

SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(SRC
    ResponseIndexer.cc
    NodeResponse.cc
    ResponseServer.cc
)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

ADD_LIBRARY(erme_response
            SHARED
            ${SRC} 
)
TARGET_LINK_LIBRARIES( erme_response
                       erme_geometry
                       serment_comm
                       ${Detran_LIBRARIES}                    
)

INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/lib erme_response)

#-------------------------------------#
# TESTING
#-------------------------------------#

if( SERMENT_ENABLE_TEST )
    include_directories(${SRC_DIR}/test)
    add_subdirectory(${SRC_DIR}/test)
endif()

