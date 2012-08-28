#------------------------------------------------------------------------------#
# serment/src/erme/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding serment erme tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_StateERME			test_StateERME.cc)
TARGET_LINK_LIBRARIES(test_StateERME 	erme)

ADD_EXECUTABLE(test_Connect 			test_Connect.cc)
TARGET_LINK_LIBRARIES(test_Connect  	erme)


#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST(test_StateERME test_StateERME 0)
