#------------------------------------------------------------------------------#
# src/pyerme/CMakeLists.txt
#------------------------------------------------------------------------------#

set_source_files_properties(serment.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(pyserment python serment.i)
SWIG_LINK_LIBRARIES(pyserment
                    ${PYTHON_LIBRARIES} 
                    serment_comm
                    erme_geometry
                    erme
                    erme_utils)

# Install the module and library
install(FILES ${CMAKE_BINARY_DIR}/src/pyerme/pyserment.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/serment)
install(FILES ${CMAKE_BINARY_DIR}/src/pyerme/_pyserment.so
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/serment)
        
# Install the helper initialization file
install_files(/python/serment .py)

# add the utilities
#add_subdirectory(pyermeutilities)

