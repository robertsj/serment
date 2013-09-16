#-------------------------------------------------------------------------------
# serment/src/xerme/CMakeLists.txt
#-------------------------------------------------------------------------------

SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(SRC
    xerme.cc
)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

add_executable(serment
               ${SRC} 
)
TARGET_LINK_LIBRARIES(serment 
                      erme_utils     
                      ${GPERFTOOLS_LIBRARIES}
)

INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/bin serment)

#-------------------------------------#
# TESTING
#-------------------------------------#

#if( SERMENT_ENABLE_TEST )
#    include_directories(${SRC_DIR}/test)
#    add_subdirectory(${SRC_DIR}/test)
#endif()

