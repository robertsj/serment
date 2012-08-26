#-------------------------------------------------------------------------------
# serment/src/comm/CMakeLists.txt
#-------------------------------------------------------------------------------

SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(SRC
    MPI.cc
    Serial.cc
)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

ADD_LIBRARY(serment_comm
            SHARED
            ${SRC} 
)
TARGET_LINK_LIBRARIES(serment_comm
                      ${Detran_LIBRARIES}
                      ${MPI_LIBRARIES}    
                      ${Boost_LIBRARIES}                 
)

INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/lib serment_comm)

#-------------------------------------#
# TESTING
#-------------------------------------#

if( SERMENT_ENABLE_TEST )
    include_directories(${SRC_DIR}/test)
    add_subdirectory(${SRC_DIR}/test)
endif()

