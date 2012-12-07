#------------------------------------------------------------------------------#
# src/pyerme/CMakeLists.txt
#------------------------------------------------------------------------------#

SET_SOURCE_FILES_PROPERTIES(serment.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(pyserment python serment.i)
SWIG_LINK_LIBRARIES(pyserment
                    ${PYTHON_LIBRARIES} 
                    serment_comm
                    erme_geometry
                    erme
                    erme_utils)

# Install the module and library
INSTALL(FILES ${CMAKE_BINARY_DIR}/src/pyerme/pyserment.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/serment)
INSTALL(FILES ${CMAKE_BINARY_DIR}/src/pyerme/_pyserment.so
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/serment)
        
# Install the helper initialization file
INSTALL_FILES(/python/serment __init__.py)

# add the utilities
#add_subdirectory(pydetranutils)

