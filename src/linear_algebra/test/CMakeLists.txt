#------------------------------------------------------------------------------#
# serment/src/linear_algebra/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding serment linear algebra tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_Vector                test_Vector.cc)
TARGET_LINK_LIBRARIES(test_Vector         linear_algebra ${Detran_LIBRARIES})

ADD_EXECUTABLE(test_Matrix                test_Matrix.cc)
TARGET_LINK_LIBRARIES(test_Matrix         linear_algebra ${Detran_LIBRARIES})

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST(test_Vector test_Vector 0)

ADD_TEST(test_Matrix test_Matrix 0)
