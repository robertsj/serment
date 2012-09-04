#------------------------------------------------------------------------------#
# serment/src/erme/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding serment erme tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_StateERME				test_StateERME.cc)
TARGET_LINK_LIBRARIES(test_StateERME 		erme)

ADD_EXECUTABLE(test_Connect 				test_Connect.cc)
TARGET_LINK_LIBRARIES(test_Connect  		erme)

ADD_EXECUTABLE(test_ResponseMatrix	    	test_ResponseMatrix.cc)
TARGET_LINK_LIBRARIES(test_ResponseMatrix 	erme)

ADD_EXECUTABLE(test_MR_scaling              test_MR_scaling.cc)
TARGET_LINK_LIBRARIES(test_MR_scaling       erme)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST(test_StateERME 		test_StateERME 0)
ADD_TEST(test_Connect   		test_Connect 0)
ADD_TEST(test_ResponseMatrix 	test_ResponseMatrix 0)
