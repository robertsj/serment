#------------------------------------------------------------------------------#
# serment/src/erme_response/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding erme_response tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_ResponseIndexer	           test_ResponseIndexer.cc)
TARGET_LINK_LIBRARIES(test_ResponseIndexer     erme_response)

ADD_EXECUTABLE(test_NodeResponse	           test_NodeResponse.cc)
TARGET_LINK_LIBRARIES(test_NodeResponse        erme_response)

ADD_EXECUTABLE(test_ResponseServer	           test_ResponseServer.cc)
TARGET_LINK_LIBRARIES(test_ResponseServer      erme_response)

ADD_EXECUTABLE(test_ResponseDatabase           test_ResponseDatabase.cc)
TARGET_LINK_LIBRARIES(test_ResponseDatabase    erme_response)

ADD_EXECUTABLE(test_Interpolation              test_Interpolation.cc)
TARGET_LINK_LIBRARIES(test_ResponseDatabase    erme_response)


#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST(test_ResponseIndexer       test_ResponseIndexer 0)
ADD_TEST(test_NodeResponse          test_NodeResponse 0)
ADD_TEST(test_ResponseServer        test_ResponseServer 0)
ADD_TEST(test_ResponseDatabase      test_ResponseDatabase 0)

#------------------------------------------------------------------------------#
# Install
#------------------------------------------------------------------------------#

message("putting test in ${CMAKE_BINARY_DIR}/src/erme_response/test ")
INSTALL(FILES test.h5 DESTINATION ${CMAKE_BINARY_DIR}/src/erme_response/test)
