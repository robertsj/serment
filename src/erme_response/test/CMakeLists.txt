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

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST(test_ResponseIndexer 	test_ResponseIndexer 0)
ADD_TEST(test_NodeResponse 		test_NodeResponse 0)
ADD_TEST(test_ResponseServer 	test_ResponseServer 0)
