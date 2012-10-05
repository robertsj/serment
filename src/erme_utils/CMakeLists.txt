#-------------------------------------------------------------------------------
# serment/src/erme_utils/CMakeLists.txt
#-------------------------------------------------------------------------------

SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(SRC
    Manager.cc
)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

ADD_LIBRARY(erme_utils
            SHARED
            ${SRC} 
)
TARGET_LINK_LIBRARIES(erme_utils
                      erme_solver                   
)

INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/lib erme_utils)

#-------------------------------------#
# TESTING
#-------------------------------------#

#if( SERMENT_ENABLE_TEST )
#    include_directories(${SRC_DIR}/test)
#    add_subdirectory(${SRC_DIR}/test)
#endif()

