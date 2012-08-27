#------------------------------------------------------------------------------#
# serment/src/comm/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding serment comm tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_PingPong                test_PingPong.cc)
TARGET_LINK_LIBRARIES(test_PingPong         serment_comm)

ADD_EXECUTABLE(test_GlobalReduction         test_GlobalReduction.cc)
TARGET_LINK_LIBRARIES(test_GlobalReduction  serment_comm)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST(test_PingPong_sendreceive 	test_PingPong 0)
ADD_TEST(test_PingPong_bandwidth 	test_PingPong 0)
ADD_TEST(test_PingPong_latency	 	test_PingPong 0)

ADD_TEST(test_GlobalReduction       test_GlobalReduction 0)
